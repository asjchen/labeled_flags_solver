# Graph Utilities

import networkx as nx 
import copy
import matplotlib.pyplot as plt

def equiv_nodes(node1, node2):
    return ((node1['flag'], node1['cluster']) == (node2['flag'], node2['cluster']))

def equiv_edges(edge1, edge2):
    return ((edge1['flags'], edge1['clusters']) == (edge2['flags'], edge2['clusters']))


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

    # TO test:
    # make sure that if we have .-.-. as a flag within a cluster, then attaching to the middle is not the same as attaching to the end
    def is_isomorphic(self, network2):
        # Both nodes and edges need to have equivalence properties
        return nx.is_isomorphic(self.graph, network2.graph, node_match=equiv_nodes,
            edge_match=equiv_edges)

    # # def check_cluster_independence?

    # def multiply(self, single_graph2):


    # def average_flags(self):

def run_tests():
    # Construct a V-graph 2-1-2, with the vertex of degree as a flag
    v_raw_graph = nx.Graph()
    v_raw_graph.add_nodes_from([1, 2, 3])
    v_raw_graph.add_edges_from([(1, 2), (1, 3)])
    cluster_dict = { 1: 1, 2: 2, 3: 2 }
    flag_dict = { 1: True, 2: False, 3: False }
    v_graph = SingleGraph(v_raw_graph, cluster_dict, flag_dict)
    v_graph.draw()

    # Reflexive Isomorphism
    print('Reflexive Isomorphism (should be true): {}'.format(v_graph.is_isomorphic(v_graph)))

    # Check isomorphism in case of edges attaching to same clusters, but fundamentally different graphs
    sim_labels_raw1 = nx.Graph()
    sim_labels_raw1.add_nodes_from([1, 2, 3, 4])
    sim_labels_raw1.add_edges_from([(1, 2), (1, 3), (3, 4)])
    sim_labels_raw2 = nx.Graph()
    sim_labels_raw2.add_nodes_from([1, 2, 3, 4])
    sim_labels_raw2.add_edges_from([(1, 2), (1, 3), (1, 4)])
    cluster_dict = { 1: 1, 2: 1, 3: 1, 4: 2 }
    flag_dict = { 1: False, 2: False, 3: False, 4: False }
    sim_labels_graph1 = SingleGraph(sim_labels_raw1, cluster_dict, flag_dict)
    sim_labels_graph2 = SingleGraph(sim_labels_raw2, cluster_dict, flag_dict)
    sim_labels_graph1.draw()
    sim_labels_graph2.draw()
    print('Same Edge Labels, Not Isomorphic (should be false): {}'.format(sim_labels_graph1.is_isomorphic(sim_labels_graph2)))

    # Check permutation of the vertices while maintaining cluster labels is isomorphic
    v_raw_graph2 = nx.Graph()
    v_raw_graph2.add_nodes_from([1, 2, 3])
    v_raw_graph2.add_edges_from([(1, 2), (2, 3)])
    cluster_dict = { 1: 2, 2: 1, 3: 2 }
    flag_dict = { 1: False, 2: True, 3: False }
    v_graph2 = SingleGraph(v_raw_graph2, cluster_dict, flag_dict)
    print('Permutation Isomorphism (should be true): {}'.format(v_graph2.is_isomorphic(v_graph)))

    # Check putting different clusters is not isomorphic
    diff_v_raw_graph = nx.Graph()
    diff_v_raw_graph.add_nodes_from([1, 2, 3])
    diff_v_raw_graph.add_edges_from([(1, 2), (1, 3)])
    cluster_dict = { 1: 1, 2: 3, 3: 3 }
    flag_dict = { 1: True, 2: False, 3: False }
    diff_v_graph = SingleGraph(diff_v_raw_graph, cluster_dict, flag_dict)
    diff_v_graph.draw()
    print('Same Graph, Different Clusters (should be false): {}'.format(diff_v_graph.is_isomorphic(v_graph)))



run_tests()

