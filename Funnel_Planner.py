import numpy as np
import matplotlib.pyplot as plt
plt.style.use("default")
from tqdm import tqdm
import pickle
# from collections import defaultdict

from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches

def make_legend_arrow(legend, orig_handle,
                      xdescent, ydescent,
                      width, height, fontsize):
    p = mpatches.FancyArrow(0, 0.5*height, width, 0, length_includes_head=True, head_width=0.75*height )
    return p

class Vertex:
    def __init__(self, position, index, funnel_id, S, rho):
        self.position = position
        self.index = index
        self.funnel_id = funnel_id
        self.S = S
        self.rho = rho
        self.next = None
#         self.neighbor_vertex = set()
        self.neighbor_vertex = []

class Funnel_Planner:
    def __init__(self, file_path="obj/funnel1"):
        with open(file_path, "rb") as f:
            self.funnel_library = pickle.load(f)
#         self.vertex_ls = defaultdict([])
        self.vertex_ls = []
        self.funnel_length = self.funnel_library[0]["re_x_values"].shape[0]
        #         self.funnel_to_vertex = []

        self.reachability_backup = {"position": np.array([np.inf, np.inf])}

    def graph_builder(self, segment_size=25, clean=True, index_range=(0, 120)):
        if clean:
            #             self.vertex_ls = defaultdict([])
            self.vertex_ls = []
#             self.funnel_to_vertex = []

        count = 0
        for i, funnel in tqdm(
                enumerate(self.funnel_library[index_range[0]:index_range[1]])):
            for j in list(range(0, self.funnel_length, segment_size)) + [-1]:
                vertex = Vertex(position=funnel["re_x_values"][j],
                                index=count,
                                funnel_id=i,
                                S=funnel["re_S_values"][j],
                                rho=funnel["rho"])
                if j != 0:
                    #                     self.vertex_ls[i][-1].next = vertex
                    #                     self.vertex_ls[i][-1].neighbor_vertex.add(vertex)
                    self.vertex_ls[-1].next = count
                    self.vertex_ls[-1].neighbor_vertex.append(count)
#                 self.vertex_ls[i].append(vertex)
                self.vertex_ls.append(vertex)
                #                 print(count, i, j)
                count += 1

        for index, vertex in enumerate(self.vertex_ls):
            vertex.neighbor_vertex += self.check_vertex_neighbor(vertex, index)

    def check_vertex_neighbor(self, vertex_, index_):
        next_vertex = []
        #         for i in range(len(self.funnel_library)):
        #             for j in range(1 + int(self.funnel_length/segment_size)):
        #                 vertex = self.vertex_ls[i][j]
        for index, vertex in enumerate(self.vertex_ls):
            if vertex.funnel_id == vertex_.funnel_id:  # don't avoid casuality
                continue
            x_bar = vertex_.position - vertex.position
            distance = x_bar.dot(vertex.S).dot(x_bar)
            if distance < vertex.rho and index != index_:
                next_vertex.append(index)

        return next_vertex

    def check_vertex_inside_final_roa(self, ):
        pass

    def check_point_inside_vertex_roa(self, position, plot=False):
        vertex_capture = []
        for index, vertex in enumerate(self.vertex_ls):
            x_bar = position - vertex.position
            distance = x_bar.dot(vertex.S).dot(x_bar)
            if distance < vertex.rho:
                vertex_capture.append(index)

        if plot:
            plt.scatter(position[0], position[1], c="g")
            for index, vertex in enumerate(self.vertex_ls):
                x_cap = vertex.position
                if index in vertex_capture:
                    plt.scatter(x_cap[0], x_cap[1], c="r")
                else:
                    plt.scatter(x_cap[0], x_cap[1], c="k")
            plt.show()

        return vertex_capture

    def check_point_reachability(self, position, iteration=2, plot=False):
        if np.isclose(self.reachability_backup["position"], position).all():
            vertex_reachable = self.reachability_backup["vertex_reachable"]
            vertex_reachable_set = self.reachability_backup[
                "vertex_reachable_set"]
            iteration_last = self.reachability_backup["iteration_last"]
        else:
            vertex_reachable = [self.check_point_inside_vertex_roa(position)]
            vertex_reachable_set = set(vertex_reachable[-1])
            iteration_last = 0

        for i in range(iteration_last, iteration):
            vertex_reachable.append([])
            for index in tqdm(vertex_reachable[i]):
                #                 vertex_capture = self.check_point_inside_vertex_roa(self.vertex_ls[index].position)
                neighbor_vertex = self.vertex_ls[index].neighbor_vertex
                for v_next in neighbor_vertex:
                    if v_next not in vertex_reachable_set:
                        vertex_reachable[-1].append(v_next)
                        while self.vertex_ls[vertex_reachable[-1][-1]].next != None \
                                and self.vertex_ls[vertex_reachable[-1][-1]].next not in vertex_reachable_set:
                            vertex_reachable[-1].append(
                                self.vertex_ls[vertex_reachable[-1][-1]].next)
#                 vertex_reachable[-1] += [self.vertex_ls[j].position for j in vertex_capture]
#                 vertex_ reachable[-1] += [v for v in vertex_capture if v not in vertex_reachable_set]
                vertex_reachable_set.update(set(vertex_reachable[-1]))
            if len(vertex_reachable[-1]) == 0:
                break
            # vertex_reachable_set.update(set(vertex_reachable[-1]))

        self.reachability_backup["position"] = position
        self.reachability_backup["vertex_reachable"] = vertex_reachable
        self.reachability_backup["vertex_reachable_set"] = vertex_reachable_set
        self.reachability_backup["iteration_last"] = iteration

        if plot:
            fig = plt.figure(figsize=(5, 5))
            plt.xlim(-2 * np.pi, 2 * np.pi)
            plt.ylim(-3.5 * np.pi, 3.5 * np.pi)
            plt.xlabel(r"$\theta$", fontsize=20)
            plt.ylabel(r"$\dot{\theta}$", fontsize=20)
            plt.scatter(position[0],
                        position[1],
                        facecolors="g",
                        edgecolors="k",
                        label="Starting Point")
            for index, vertex in enumerate(self.vertex_ls):
                x_cap = vertex.position
                if index in vertex_reachable_set:
                    plt.scatter(x_cap[0],
                                x_cap[1],
                                facecolors="r",
                                edgecolors="k")
                else:
                    plt.scatter(x_cap[0],
                                x_cap[1],
                                facecolors="None",
                                edgecolors="k")

            # dummy
            plt.scatter(100,
                        100,
                        facecolors="r",
                        edgecolors="k",
                        label="Reachable Vertices in {}st iteration".format(iteration+1))
            plt.scatter(100,
                        100,
                        facecolors="None",
                        edgecolors="k",
                        label="Nonreachable Vertices in {}st iteration".format(iteration+1))

            plt.legend(loc="best")
            # plt.savefig("fig/reachability_{}.png".format(iteration),
            #             dpi=300,
            #             bbox_inches='tight')
            
            plt.show()

        return vertex_reachable_set

    def graph_visualizer(self,
                         x0lim=np.pi,
                         x1lim=np.pi * 3.5,
                         figsize=8,
                         no_interconnect=False):
        edge_funnel = []
        edge_inter = []

        fig = plt.figure(figsize=(5, 5))
        #             fig = plt.figure(figsize=(figsize, figsize * x1lim / x0lim))
        plt.xlim(-x0lim, x0lim)
        plt.ylim(-x1lim, x1lim)
        plt.xlabel(r"$\theta$", fontsize=20)
        plt.ylabel(r"$\dot{\theta}$", fontsize=20)

        for index, vertex in enumerate(self.vertex_ls):
            curr_x, curr_y = vertex.position
            plt.scatter(curr_x, curr_y, facecolors='none', edgecolors='k')
            if vertex.next != None:
                edge_funnel.append((index, vertex.next))
            for neighbor_index in vertex.neighbor_vertex[1:]:
                edge_inter.append((index, neighbor_index))

        for e in tqdm(edge_funnel):
            if np.isclose(self.vertex_ls[e[1]].position,
                          self.vertex_ls[e[0]].position).all():
                continue
            plt.arrow(*self.vertex_ls[e[0]].position,
                      *(self.vertex_ls[e[1]].position -
                        self.vertex_ls[e[0]].position),
                      head_width=0.2,
                      head_length=0.4,
                      fc='r',
                      ec='b',
                      shape="left",
                      length_includes_head=True)

        if not no_interconnect:
            for e in tqdm(edge_inter):
                if np.isclose(self.vertex_ls[e[1]].position,
                              self.vertex_ls[e[0]].position).all():
                    continue
                plt.arrow(*self.vertex_ls[e[0]].position,
                          *(self.vertex_ls[e[1]].position -
                            self.vertex_ls[e[0]].position),
                          head_width=0.2,
                          head_length=0.4,
                          fc='r',
                          ec='g',
                          shape="left",
                          length_includes_head=True)
        # dummy
        scatter = plt.scatter(100,
                    100,
                    facecolors='none',
                    edgecolors='k',
                    label="Vertices")
        arrow1 = plt.arrow(100,
                           100,
                           0,
                           0,
                           head_width=0.2,
                           head_length=0.4,
                           fc='r',
                           ec='b',
                           shape="left",
                           length_includes_head=True,
                           label="Edges inside Funnels")
        if not no_interconnect:
            arrow2 = plt.arrow(100,
                               100,
                               0,
                               0,
                               head_width=0.2,
                               head_length=0.4,
                               fc='r',
                               ec='g',
                               shape="left",
                               length_includes_head=True,
                               label="Interconnective Edges")
            plt.legend([scatter,arrow1,
                        arrow2,], 
                        ['Vertices','Edges inside Funnels',
                        'Interconnective Edges',],
                        handler_map={
                           mpatches.FancyArrow:
                           HandlerPatch(patch_func=make_legend_arrow),
                        },
                       loc="best"
            )
        else:
            plt.legend([scatter,arrow1], 
                        ['Vertices','Edges inside Funnels'],
                        handler_map={
                           mpatches.FancyArrow:
                           HandlerPatch(patch_func=make_legend_arrow),
                        },
                       loc="best"
            )
#         plt.legend(loc="best")

        plt.show()
        # if no_interconnect:
        #     figure_name = "graph_no_interconnect"
        # else:
        #     figure_name = "graph_w_interconnect"
        # plt.savefig("fig/{}.png".format(figure_name), dpi=300, bbox_inches='tight')


#     def funnel_visualizer_2d(self,):
#             fig = plt.figure(figsize=(5,5))
#             plt.xlim(-2*np.pi, 2*np.pi)
#             plt.ylim(-4*np.pi, 4*np.pi)