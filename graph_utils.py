# Graph Utilities

import networkx as nx 
import copy
import matplotlib.pyplot as plt

class SingleGraph:
    def __init__(self, orig_graph, cluster_dict, flag_dict):
        self.graph = copy.deepcopy(orig_graph)

        # Check cluster and flag labels are correct, and then add them
        try:
            assert list(self.graph.nodes) == cluster_dict.keys()
        except AssertionError as err:
            print('Cluster labels ({}) inconsistent with provided graph.'.format(cluster_dict.keys()))
            raise 
        try:
            assert list(self.graph.nodes) == flag_dict.keys()
        except AssertionError as err:
            print('Flag labels ({}) inconsistent with provided graph.'.format(flag_dict.keys()))
            raise 
        for node_idx in cluster_dict:
            self.graph.add_node(node_idx, cluster=cluster_dict[node_idx],
                flag=flag_dict[node_idx])
        
        # Add edge and flag attributes as sets of cluster labels
        for n1, n2 in self.graph.edges:
            self.graph.add_edge(n1, n2, 
                clusters=set([cluster_dict[n1], cluster_dict[n2]]),
                flags=set([flag_dict[n1], flag_dict[n2]]))

    def draw(self, flag_size=400, node_size=200):
        pos = nx.spectral_layout(self.graph)
        nx.draw_networkx_nodes(self.graph, pos, node_size=node_size)

        flag_status = nx.get_node_attributes(self.graph, 'flag')
        flags = [n for n in flag_status if flag_status[n]]
        nx.draw_networkx_nodes(self.graph, pos, nodelist=flags, 
            node_shape='d', node_size=flag_size)

        cluster_labels = nx.get_node_attributes(self.graph, 'cluster')
        nx.draw_networkx_labels(self.graph, pos, labels=cluster_labels)

        nx.draw_networkx_edges(self.graph, pos)
        plt.show()
        plt.close()


    # def equiv_nodes(self, node1, node2):
    #     return ((node1['flag'], node1['cluster']) == (node2['flag'], node2['cluster']))

    # def equiv_edges(self, edge1, edge2):
    #     return ((edge1['flags'], edge1['clusters']) == (edge2['flags'], edge2['clusters']))

    # # TO test:
    # # make sure that if we have .-.-. as a flag within a cluster, then attaching to the middle is not the same as attaching to the end
    # def is_isomorphic(self, network2):
    #     # Both nodes and edges need to have equivalence properties
    #     return nx.is_isomorphic(self.graph, network2.graph, node_match=equiv_nodes,
    #         edge_match=equiv_edges)

    # # def check_cluster_independence?

    # def multiply(self, single_graph2):


    # def average_flags(self):

def run_tests():
    # Construct a V-graph 2-1-2, with the vertex of degree as a flag
    v_raw_graph = nx.Graph()
    v_raw_graph.add_node(1)
    v_raw_graph.add_node(2)
    v_raw_graph.add_node(3)
    v_raw_graph.add_edge(1, 2)
    v_raw_graph.add_edge(1, 3)
    cluster_dict = { 1: 1, 2: 2, 3: 2 }
    flag_dict = { 1: True, 2: False, 3: False }
    v_graph = SingleGraph(v_raw_graph, cluster_dict, flag_dict)
    v_graph.draw()

run_tests()

