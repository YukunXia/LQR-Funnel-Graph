from Funnel_Planner import *

planner = Funnel_Planner(file_path="obj/funnel1.obj")
planner.graph_builder(segment_size=15)

vertex_reachable_set = planner.check_point_reachability(position=np.array([0,0]), iteration=5, plot=True)