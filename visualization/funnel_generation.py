import sys
sys.path.append("..")

from Pendulum import *
LOAD = True

pendulum = Pendulum()

if not LOAD:
    pendulum.max_dt = 0.2
    pendulum.R = np.array([[2]])

    pendulum.funnel_library_generator(total_funnel=120, least_total_steps=5, LQR=True, save=False)

    pendulum.R = np.array([[20]])
    pendulum.max_dt = 0.5
else:
    pendulum.funnel_library = pickle.load(open("obj/funnel1.obj", "rb"))

pendulum.funnel_library_visualizer(start_index=10, end_index=30)