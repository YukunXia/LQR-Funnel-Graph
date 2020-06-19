from pydrake.all import (DirectCollocation, PiecewisePolynomial, Solve, eq,
                         FirstOrderTaylorApproximation,
                         LinearQuadraticRegulator,
                         RealContinuousLyapunovEquation, MathematicalProgram,
                         Variables, Solve, sin)  # , Linearize)
from pydrake.examples.pendulum import PendulumPlant, PendulumState, PendulumParams

from underactuated.pendulum import PendulumVisualizer
from underactuated import plot_2d_phase_portrait

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import mpl_toolkits.mplot3d.art3d as art3d
plt.style.use("default")

from tqdm import tqdm
import pickle

class Pendulum:
    def __init__(self, ):
        self.pendulum_plant = PendulumPlant()
        self.pendulum_context = self.pendulum_plant.CreateDefaultContext()
        self.mass, self.length, self.damping, self.gravity = PendulumParams(
        ).CopyToVector()

        self.R = np.array([[20]])
        self.R_inv = np.linalg.inv(self.R)
        self.Q = np.diag([10, 1])

        self.min_dt = 0.05
        self.max_dt = 0.5
        self.torque_limit = 3.0

        self.solve_lyapunov_eq = False

        self.x_trajectory = None
        self.u_trajectory = None
        self.times = None
        self.x_values = None
        self.u_lookup = None
        self.u_values = None
        self.S_values = None
        self.truncate_index = None
        self.funnel_library = None
        self.common_rho = None

        self.vertex_ls = None

    def direct_collocation(self,
                           x_init=np.array([0.0, 0.0]),
                           x_final=np.array([np.pi, 0.0]),
                           N=21):
        #         N: dircol breaks
        #         min_dt: 0.05
        #         max_dt: 0.5
        #         dt: delta t
        #         tf: time final \in [N*min_dt, N*max_dt] (N-1 here?)

        dircol = DirectCollocation(self.pendulum_plant,
                                   self.pendulum_context,
                                   num_time_samples=N,
                                   minimum_timestep=self.min_dt,
                                   maximum_timestep=self.max_dt)
        # 1. Constraints

        dircol.AddEqualTimeIntervalsConstraints()

        u = dircol.input()  # size depends on the plant
        dircol.AddConstraintToAllKnotPoints(-self.torque_limit <= u[0])
        dircol.AddConstraintToAllKnotPoints(u[0] <= self.torque_limit)

        initial_state = PendulumState()
        initial_state.set_theta(x_init[0])
        initial_state.set_thetadot(x_init[1])

        dircol.AddLinearConstraint(
            eq(dircol.initial_state(), initial_state.get_value()))

        final_state = PendulumState()
        final_state.set_theta(x_final[0])
        final_state.set_thetadot(x_final[1])

        dircol.AddLinearConstraint(
            eq(dircol.final_state(), final_state.get_value())
        )  # state is a vector, use eq() rather than ==

        # 2. Cost

        dircol.AddRunningCost(self.R[0, 0] * u[0]**2)

        dircol.AddFinalCost(dircol.time())

        # 3. Init traj

        initial_x_trajectory = PiecewisePolynomial.FirstOrderHold(
            [0., 4.],
            [initial_state.get_value(),
             final_state.get_value()])  # first arg is the break points
        # i.e. t = 0 -> x = initial_state; t = 4 -> x = final_state
        # https://drake.mit.edu/doxygen_cxx/classdrake_1_1trajectories_1_1_piecewise_polynomial.html
        # A note on terminology. For piecewise-polynomial interpolation,
        # we use breaks to indicate the scalar (e.g. times) which form the boundary of each segment.
        # We use samples to indicate the function value at the breaks, e.g. p(breaks[i]) = samples[i].
        # The term knot should be reserved for the "(x,y)" coordinate, here knot[i] = (breaks[i], samples[i]),
        # though it is used inconsistently in the interpolation literature (sometimes for breaks,
        # sometimes for samples), so we try to mostly avoid it here.

        dircol.SetInitialTrajectory(PiecewisePolynomial(),
                                    initial_x_trajectory)
        # 4. Solve
#         print("\nbefore dircol solve")
#         print("x_init =", x_init)
#         print("x_final =", x_final)
        result = Solve(dircol)
        #         assert result.is_success()
#         print("result =", result.is_success())
        # 5. Reconstruct continuous x and u

        self.x_trajectory = dircol.ReconstructStateTrajectory(result)
        self.u_trajectory = dircol.ReconstructInputTrajectory(result)

        return result.is_success()

    def traj_discrete(self, size=100):
        self.times = np.linspace(self.x_trajectory.start_time(),
                                 self.x_trajectory.end_time(), size)
        self.x_values = np.vstack(
            [self.x_trajectory.value(t).T for t in self.times])

        self.u_lookup = np.vectorize(self.u_trajectory.value)
        self.u_values = self.u_lookup(self.times)

    def linearization(self, x, u):
        local_context = self.pendulum_plant.CreateDefaultContext()

        local_context.get_mutable_continuous_state_vector().\
            SetFromVector(x)
        self.pendulum_plant.get_input_port().FixValue(local_context, u)
        input_i = self.pendulum_plant.get_input_port().get_index(
        )  # input_i = 3?
        output_i = self.pendulum_plant.get_state_output_port().get_index(
        )  # output_i = 1?
        pendulum_lin = FirstOrderTaylorApproximation(
            self.pendulum_plant,
            local_context,
            input_port_index=input_i,
            output_port_index=output_i)
        return pendulum_lin.A(), pendulum_lin.B()

    def LQR_get_K_P(self, x, u):
        A, B = self.linearization(x, u)
        K, S = LinearQuadraticRegulator(A, B, self.Q, self.R)
        if self.solve_lyapunov_eq:
            A_new = A - B @ K
            P = RealContinuousLyapunovEquation(A_new, np.eye(2))
        else:
            P = S
        return K, P

    def f_local(self, x_bar, x_0, K, symoblic=False):
        q_bar, q_dot_bar = x_bar

        u_bar = -K[0][0] * q_bar - K[0][1] * q_dot_bar
        #         u_bar = - K.dot(x_bar)
        # u_bar = -K * x_bar, u = u_bar + 0
        if symoblic:
            x_dot_bar = [
                q_dot_bar,
                (
                    u_bar - self.damping * q_dot_bar -
                    self.mass * self.gravity * self.length *
                    (np.cos(x_0[0]) * q_bar - (np.sin(x_0[0]) * q_bar**2) / 2 -
                     (np.cos(x_0[0]) * q_bar**3) / 6) +
                    (np.sin(x_0[0]) * q_bar**4) / 24
                    #                      - (np.cos(x_0[0]) * q_bar**5) / 120
                ) / (self.mass * self.length**2)
                # 4th order taylor expansion of sin
            ]
        else:
            x_dot_bar = np.array([
                q_dot_bar,
                (u_bar - self.damping * q_dot_bar -
                 self.mass * self.gravity * self.length *
                 (np.sin(q_bar + x_0[0])) - np.sin(x_0[0])) /
                (self.mass * self.length**2)
            ])
        return x_dot_bar


    def traj_generator(self,
                       num=3,
                       x0lim=np.pi,
                       x1lim=np.pi * 3.5,
                       growth_len=np.pi / 6,
                       plot=True,
                       size=100,
                       deterministic=False,
                       x_init_given=None,
                       x_final_given=None
                       ):
        if plot:
            fig = plt.figure(figsize=(5,5))
            ax = fig.gca()
            ax.scatter(20,
                       20,
                       facecolors="r",
                       edgecolors="none",
                       label="Starting Point")
            ax.scatter(20,
                       20,
                       facecolors="k",
                       edgecolors="none",
                       label="Ending Point")
        iterator = tqdm(range(num)) if plot else range(num)
        if not deterministic:
            for i in iterator:
                result = False
                while result == False:
                    x_init = (np.random.rand((2)) - 0.5) * 2
                    x_init[0] *= x0lim
                    x_init[1] *= x1lim

                    theta = np.random.uniform(0, 2 * np.pi)
                    x_final = x_init + growth_len * np.array(
                        [np.sin(theta), np.cos(theta)])

                    result = self.direct_collocation(x_init, x_final)
                    self.traj_discrete(size)
                if plot:
                    linewidth = 0.1
                    ax.plot(self.x_values[:, 0],
                            self.x_values[:, 1],
                            linewidth=linewidth)
                    ax.scatter(self.x_values[0, 0],
                               self.x_values[0, 1],
                               facecolors="r",
                               edgecolors="none")
                    ax.scatter(self.x_values[-1, 0],
                               self.x_values[-1, 1],
                               facecolors="k",
                               edgecolors="none")
        else:
            result = self.direct_collocation(x_init_given, x_final_given)
#             print("finish dir col in traj generation in deterministic mode")
            if not result:
                return False
            self.traj_discrete(size)
            if plot:
                ax.plot(self.x_values[:, 0],
                        self.x_values[:, 1],
                        linewidth=0.1)
                ax.scatter(self.x_values[0, 0],
                           self.x_values[0, 1],
                           facecolors="r",
                           edgecolors="none")
                ax.scatter(self.x_values[-1, 0],
                           self.x_values[-1, 1],
                           facecolors="k",
                           edgecolors="none")
        if plot:
            ax.set_xlabel(r"$\theta$", fontsize=20)
            ax.set_ylabel(r"$\dot{\theta}$", fontsize=20)
            ax.set_xlim(-x0lim*2, x0lim*2)
            ax.set_ylim(-x1lim, x1lim)
            ax.legend(fontsize=10, loc="best")
            # plt.savefig("random_trajectories.png", dpi=300, bbox_inches = 'tight')
            plt.show()

        return True

    def traj_rediscrete(self, total_time, size=100, LQR=True):
        end_time = self.x_trajectory.end_time(
        ) if LQR else self.x_trajectory.end_time() * (size - 1) / size
        start_time = end_time - total_time
        re_times = np.linspace(start_time, end_time, size)
        re_x_values = np.vstack(
            [self.x_trajectory.value(t).T for t in re_times])
        re_u_values = self.u_lookup(re_times)
        
        re_K_values = []
        for i in range(size):
            re_K_values.append(self.LQR_get_K_P(re_x_values[i], re_u_values[i])[0])
        re_K_values = np.array(re_K_values)
    
        t_span = [0, total_time]
        t_eval = np.linspace(*t_span, size)

        P = self.S_values[-1] if LQR else self.S_values[-2]

        sol = solve_ivp(fun=self.differential_riccati_eq,
                        t_span=t_span,
                        y0=P.flatten(),
                        method="RK45",
                        t_eval=t_eval)
        re_S_values = sol.y.T[::-1, :].reshape(-1, 2, 2)
        for S in re_S_values:
            assert self.is_pos_def(S) == True

        return re_x_values, re_u_values, re_S_values, re_K_values

    # good means x_f is outside of the terminal ellipses of other funnels
    def is_a_good_traj(self, ):
        current_end = self.x_values[-1]
        for funnel in self.funnel_library:
            funnel_end = funnel["x_values"][-1]
            funnel_rho = funnel["rho"]
            funnel_S = funnel["S_values"][-1]
            x_bar = current_end - funnel_end
            x_bar_S_x_bar = x_bar.dot(funnel_S).dot(x_bar)
            if x_bar_S_x_bar < funnel_rho:
                return False
        return True

    def differential_riccati_eq(self, t, S_vec):
        #         print(t, S_vec)
        length = int(np.sqrt(S_vec.size))
        S_matrix = S_vec.reshape(length, length)

        x_curr = self.x_trajectory.value(self.x_trajectory.end_time() -
                                         t)  # (n, 1)
        u_lookup = np.vectorize(self.u_trajectory.value)
        u_curr = u_lookup(self.x_trajectory.end_time() - t)  # (1, 1) -> float
        A, B = self.linearization(x_curr, u_curr)

        dS_negative = self.Q - S_matrix @ B @ self.R_inv @ B.T @ S_matrix + \
            S_matrix @ A + A.T @ S_matrix
        return dS_negative.flatten()

    def is_pos_def(self, M):
        return np.all(np.linalg.eigvals(M) > 0)

    def solve_S_traj(self, size=100, LQR=True, P=np.diag([3, 3])):
        assert self.x_trajectory.start_time()**2 < 0.0001
        t_span = [0, self.x_trajectory.end_time()]
        t_eval = np.linspace(*t_span, size)
        if LQR:  # otherwise, P is given by argument
            _, P = self.LQR_get_K_P(self.x_values[-1], self.u_values[-1])

        sol = solve_ivp(fun=self.differential_riccati_eq,
                        t_span=t_span,
                        y0=P.flatten(),
                        method="RK45",
                        t_eval=t_eval)
        self.S_values = sol.y.T[::-1, :].reshape(-1, 2, 2)
        for S in self.S_values:
            assert self.is_pos_def(S) == True

    def S_procedure_fixed_point(self,
                                x_0,
                                u_0,
                                P_given_flag=False,
                                P_given=None):
        K, P = self.LQR_get_K_P(x_0, u_0)
        if P_given_flag:
            P = P_given
        # initialize optimization problem
        prog = MathematicalProgram()
        # SOS indeterminates
        x_bar = prog.NewIndeterminates(2, 'x_bar')
        # Lyapunov function
        V = x_bar.dot(P).dot(x_bar)
        V_dot = 2 * x_bar.dot(P).dot(self.f_local(x_bar, x_0, K,
                                                  symoblic=True))
        # degree of the polynomial lambda(x)
        # no need to change it, but if you really want to,
        # keep l_deg even and do not set l_deg greater than 10
        l_deg = 4
        assert l_deg % 2 == 0
        # SOS Lagrange multipliers
        if P_given_flag:
            l = prog.NewSosPolynomial(Variables(x_bar),
                                      l_deg)[0].ToExpression()
        else:
            l = prog.NewFreePolynomial(Variables(x_bar), l_deg).ToExpression()
        # level set as optimization variable
        rho = prog.NewContinuousVariables(1, 'rho')[0]
        # write here the SOS condition described in the "Not quite there yet..." section above
        prog.AddSosConstraint(x_bar.dot(x_bar) * (V - rho) - l * V_dot)
        # insert here the objective function (maximize rho)
        prog.AddLinearCost(-rho)
        # solve program only if the lines above are filled
        if len(prog.GetAllConstraints()) != 0:
            # solve SOS program
            result = Solve(prog)
            # get maximum rho
            assert result.is_success()
            rho = result.GetSolution(rho)
            # print maximum rho

#             print(f'Verified rho = {rho}.')
        return rho

    def S_procedure_traj(self,
                         S_mode="piecewise_const",
                         rho_mode="const",
                         terminal_rho=None,
                         LQR=True):
        # run this function after running pendulum.solve_S_traj
        size = len(self.x_values)
        # time_interval = self.x_trajectory.end_time() / size

        if S_mode == "piecewise_const" and rho_mode == "const":
            rho_common = np.inf
            start_index = -1 if LQR else -2
            for i in range(start_index, -size - 1, -1):
                rho = self.S_procedure_fixed_point(self.x_values[i],
                                                   self.u_values[i],
                                                   True,
                                                   P_given=self.S_values[i])
                #                 print(rho)
                if rho < 0.5:
                    self.truncate_index = i + 1
                    return rho_common
                    break
                rho_common = min(rho_common, rho)
#                 print(rho, rho_common)
            self.truncate_index = -size
            return rho_common
        else:
            raise NotImplementedError("unimplemented mode")


    def funnel_library_generator(self,
                                 total_funnel=5,
                                 least_total_steps=5,
                                 LQR=True,
                                 save=False):  # 0.1*5
        count = 0
#         num_back_to_non_deterministic = 0
        deterministic_successful = 0
        deterministic_fail_on_dircol = 0
        deterministic_fail_on_short = 0
        self.funnel_library = []
        with tqdm(total=total_funnel) as pbar:
            deterministic = False
            x_init_given = None
            x_final_given = None
            while count < total_funnel:
                result = self.traj_generator(num=1,
                                             growth_len=np.pi,
                                             deterministic=deterministic,
                                             x_init_given=x_init_given,
                                             x_final_given=x_final_given)
#                 print("mode is deterministic =", deterministic,", result is", result)
                if deterministic and not result:
                    deterministic_fail_on_dircol += 1
                if ( not result ) or ( count > int(total_funnel / 3) and not self.is_a_good_traj() ):
                    if deterministic:
                        deterministic = False
                        x_init_given = None
                        x_final_given = None
                    continue
#                 if not self.is_a_good_traj():
#                     continue
                self.solve_S_traj(size=100, LQR=LQR)
                self.common_rho = self.S_procedure_traj(
                    S_mode="piecewise_const", rho_mode="const", LQR=LQR)
                total_steps = (-int(LQR == False) - self.truncate_index)
                total_time = total_steps * self.x_trajectory.end_time() / 100 
                if total_steps >= least_total_steps:
                    if deterministic:
                        deterministic_successful += 1
                    re_x_values, re_u_values, re_S_values, re_K_values = self.traj_rediscrete(
                        total_time, size=100, LQR=LQR)
                    self.funnel_library.append({
                        "S_values":
                        self.S_values[100 + self.truncate_index:100],
                        "x_values":
                        self.x_values[100 + self.truncate_index:100],
                        "u_values":
                        self.u_values[100 + self.truncate_index:100],
                        "rho":
                        self.common_rho,
                        #                         "total_time":
                        #                         total_time,
                        "time_interval":
                        self.x_trajectory.end_time() / 100,
                        "total_steps":
                        -self.truncate_index,
                        "re_x_values":
                        re_x_values,
                        "re_u_values":
                        re_u_values,
                        "re_S_values":
                        re_S_values,
                        "re_K_values":
                        re_K_values
                    })
                    count += 1
                    pbar.update(1)
                elif deterministic:
                    deterministic_fail_on_short += 1

                # redo the direct collocation
                if self.truncate_index > -least_total_steps or self.common_rho == np.inf or \
                    self.truncate_index <= -self.x_values.shape[0] + least_total_steps:
                    deterministic = False
                    x_init_given = None
                    x_final_given = None
                else:
                    deterministic = True
                    x_init_given = self.x_values[0]
                    x_final_given = self.x_values[self.truncate_index]

        pbar.close()
        print("deterministic ratio =", deterministic_successful/total_funnel)
        print("deterministic fail on dir col / success", deterministic_fail_on_dircol/deterministic_successful)
        print("deterministic fail on short / success", deterministic_fail_on_short/deterministic_successful)
        
        if save:
            pickle.dump(self.funnel_library, open("obj/funnel1.obj", "wb"))

    def funnel_library_visualizer(self,
                                  rebuilt=True,
                                  start_index=0,
                                  end_index=-1):
        def orientation_from_covariance(cov, sigma):
            vals, vecs = np.linalg.eig(cov)
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            w, h = 2 * np.sqrt(rho / vals)
            return w, h, theta

        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel(r"$\theta$", fontsize=10)
        ax.set_ylabel(r"$\dot{\theta}$", fontsize=10)
        ax.set_zlabel(r"$t'$", fontsize=10)

        c1 = np.array([1, 0, 0])
        c2 = np.array([0, 0, 1])

        for funnel_discription in tqdm(self.funnel_library[
                start_index:end_index]):
            rho = funnel_discription["rho"]
            if rebuilt:
                total_steps = 100
            else:
                total_steps = funnel_discription["total_steps"]

            for i in range(total_steps):
                if rebuilt:
                    cov = funnel_discription["re_S_values"][i]
                    mu = funnel_discription["re_x_values"][i]
                else:
                    cov = funnel_discription["S_values"][i]
                    mu = funnel_discription["x_values"][i]

                x, y, angle = orientation_from_covariance(cov, rho)
                e = Ellipse(mu, x, y, angle=angle)
                e.set_alpha(0.5)
                e.set_linewidth(2)
                color = i / total_steps * c1 + (total_steps -
                                                i) / total_steps * c2
                e.set_edgecolor(color)
                e.set_facecolor(color)
                e.set_fill(False)

                ax.add_patch(e)

                art3d.pathpatch_2d_to_3d(e,
                                         z=(total_steps - 1 - i) / total_steps,
                                         zdir="z")

        ax.set_xlim(-5, 5)
        ax.set_ylim(-10, 10)
        ax.set_zlim(0, 1)
        # plt.savefig("fig/funnel120.png", dpi=300, bbox_inches = 'tight')

        plt.show()
            

    # direct conection starts from here

    class Vertex:
        def __init__(self, position):
            self.position = position
            self.neighbors = set()
            self.cost_to_neighbors = dict()

    def graph_builder(
        self,
        mode="given",
        num=120*8,
        x0lim=np.pi,
        x1lim=np.pi * 3.5,
        #                       growth_len=np.pi / 6,
        plot=False):
        if mode == "given":
            np.random.seed(2020)
            points = (np.random.rand(num, 2)-0.5)*2
            np.random.seed()
            points[:, 0] *= x0lim
            points[:, 1] *= x1lim

            tri = Delaunay(points)

            edge_set = set()
            for simplex in tri.simplices:
                sorted_simplex = sorted(simplex)
                edge_set.add((sorted_simplex[0], sorted_simplex[1]))
                edge_set.add((sorted_simplex[0], sorted_simplex[2]))
                edge_set.add((sorted_simplex[1], sorted_simplex[2]))

            self.vertex_ls = [self.Vertex(points[i]) for i in range(num)]
            for simplex in tri.simplices:
                for vertex in simplex:
                    other = set(simplex)
                    other.remove(vertex)
                    self.vertex_ls[vertex].neighbors = self.vertex_ls[
                        vertex].neighbors.union(other)
            if plot:
                fig = plt.figure(figsize=(8, 16 * x1lim / x0lim))
                plt.xlim(0, x0lim)
                plt.ylim(0, x1lim)

            for edge in tqdm(edge_set):
                vertex0 = edge[0]
                vertex1 = edge[1]
                successful = self.direct_collocation(x_init=points[vertex0],
                                                     x_final=points[vertex1])
                if successful:
                    self.traj_discrete()
                    cost = self.x_trajectory.end_time() * np.average(
                        self.u_values**2) * self.R[0, 0]
                    self.vertex_ls[vertex0].cost_to_neighbors[vertex1] = cost
                    if plot:
                        plt.arrow(*points[vertex0],
                                  *(points[vertex1] - points[vertex0]),
                                  head_length=0.05,
                                  fc='g',
                                  ec='g',
                                  shape="left",
                                  length_includes_head=True)

            if plot:
                self.graph_visualizer(mode, x0lim, x1lim)

        elif mode == "incremental":
            pass
        else:
            raise NotImplementedError(
                "graph building mode {} undefined".format(mode))

    def graph_visualizer(self,
                         mode="given",
                         x0lim=np.pi,
                         x1lim=np.pi * 3.5,
                         figsize=5):
        if mode == "given":
            single_edge = set()
            double_edge = set()

            fig = plt.figure(figsize=(figsize, figsize))
            plt.xlim(-x0lim, x0lim)
            plt.ylim(-x1lim, x1lim)
            plt.xlabel(r"$\theta$", fontsize=20)
            plt.ylabel(r"$\dot{\theta}$", fontsize=20)

            for curr, vertex in enumerate(self.vertex_ls):
                curr_x, curr_y = vertex.position
                plt.scatter(curr_x, curr_y, facecolors='none', edgecolors='k')
                for neighbor in vertex.cost_to_neighbors:
                    if (neighbor, curr) in single_edge:
                        single_edge.remove((neighbor, curr))
                        double_edge.add((min(curr,
                                             neighbor), max(curr, neighbor)))
                    else:
                        single_edge.add((curr, neighbor))

            print("single edge size = {}".format(len(single_edge)))
            for edge in single_edge:
                plt.arrow(*self.vertex_ls[edge[0]].position,
                          *(self.vertex_ls[edge[1]].position -
                            self.vertex_ls[edge[0]].position),
                          head_width=0.05,
                          head_length=0.1,
                          fc='r',
                          ec='b',
                          shape="left",
                          length_includes_head=True)

            print("double edge size = {}".format(len(double_edge)))
            for edge in double_edge:
                p0x, p0y = self.vertex_ls[edge[0]].position
                p1x, p1y = self.vertex_ls[edge[1]].position
                plt.plot([p0x, p1x], [p0y, p1y], c="r")
            plt.show()

        else:
            raise NotImplementedError(
                "graph visualzation mode {} undefined".format(mode))
            
        # plt.savefig("fig/direct_connect.png", dpi=300, bbox_inches = 'tight')