import sys
sys.path.append("..")

from Funnel_Planner import *

NO_INTERCONNECT = False

planner = Funnel_Planner(file_path="obj/funnel1.obj")
planner.graph_builder(segment_size=15)

planner.graph_visualizer(figsize=10, no_interconnect=NO_INTERCONNECT)