from Pendulum import *

pendulum = Pendulum()

# pendulum.R = np.array([[0.1]])
# pendulum.max_dt = 0.1

pendulum.direct_collocation()
pendulum.traj_discrete(size = 100)

fig, [ax0, ax1] = plt.subplots(figsize=(10, 5), nrows=1, ncols=2)

ax0.plot(pendulum.x_values[:, 0], pendulum.x_values[:, 1])
ax0.set_ylabel(r"$\dot{\theta}$", fontsize=20)
ax0.set_xlabel(r"$\theta$", fontsize=20)

ax1.plot(pendulum.times, pendulum.u_values)
ax1.set_ylabel(r"Input $u(t)$", fontsize=20)
ax1.set_xlabel(r"Time $t$", fontsize=20)

plt.show()


# run the following code if in jupyter notebook

# from IPython.display import HTML

# # Animate the result.
# vis = PendulumVisualizer(show=False)
# ani = vis.animate(pendulum.x_trajectory)
# HTML(ani.to_jshtml())