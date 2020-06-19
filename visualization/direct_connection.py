import sys
sys.path.append("..")

from Pendulum import *

pendulum = Pendulum()
pendulum.graph_builder(num=120*8) # compared with 120 funnels and 8 points per funnel
pendulum.graph_visualizer()