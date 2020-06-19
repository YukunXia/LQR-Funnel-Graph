import sys
sys.path.append("..")

from Pendulum import *

pendulum = Pendulum()

pendulum.traj_generator(num=200, growth_len=np.pi)