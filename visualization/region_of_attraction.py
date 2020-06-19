import sys
sys.path.append("..")

from Pendulum import *

pendulum = Pendulum()

# get a fixed point
# theta = np.pi * 3 / 3
theta = np.pi
u = pendulum.mass * pendulum.gravity * pendulum.length * np.sin(theta)

# use the LQR solve S as the P of Lyapunov function
pendulum.solve_lyapunov_eq = False

x_0=np.array([theta, 0])
u_0=u

rho = pendulum.S_procedure_fixed_point(x_0, u_0)

# given rho plots the boundary
# of the the set L(rho) defined above
x1lim=(-5, 5, 100)
x2lim=(-20, 20, 100)
plot_torque_limit=True

K, P = pendulum.LQR_get_K_P(x_0, u_0)

# grid of the state space
x1 = np.linspace(*x1lim)
x2 = np.linspace(*x2lim)
X1, X2 = np.meshgrid(x1, x2)

# function that evaluates V(x) at a given x
# (looks bad, but it must accept meshgrids)
def eval_V(x):
    return sum(
        sum(x[i] * x[j] * Pij for j, Pij in enumerate(Pi))
        for i, Pi in enumerate(P))

# contour plot with only the rho level set
plt.contour(X1,
            X2,
            eval_V([X1, X2]),
            levels=[rho],
            colors='r',
            linewidths=3,
            zorder=3)

# fake plot for legend
plt.plot(0,
        0,
        color='r',
        linewidth=3,
        label=r'$\{ \mathbf{x} : V(\mathbf{x}) = \rho \}$')

# function that evaluates Vdot(x) at a given x
def eval_Vdot(x):
    return 2 * sum(
        sum(x[i] * pendulum.f_local(x, x_0, K)[j] * Pij
            for j, Pij in enumerate(Pi)) for i, Pi in enumerate(P))

# contour plot with only the rho level set
cs1 = plt.contour(X1,
                    X2,
                    eval_Vdot([X1, X2]),
                    colors='b',
                    levels=np.linspace(-1000, 1000, 9))
plt.gca().clabel(cs1, inline=1, fontsize=10)

# misc plot settings
plt.xlabel(r'$x_1$', fontsize=20)
plt.ylabel(r'$x_2$', fontsize=20)

# fake plot for legend
plt.plot(0, 0, color='b', label=r'$\dot{V}(\mathbf{x})$')

# plot torque limit
if plot_torque_limit:
    plt.plot(x1, (pendulum.torque_limit + K[0][0] * x1) / (-K[0][1]),
                c="g",
                label="torque limit")
    plt.plot(x1, (-pendulum.torque_limit + K[0][0] * x1) / (-K[0][1]),
                c="g")

plt.legend(loc="best")

plt.show()
