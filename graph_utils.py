# Graph Utilities

import math
import networkx as nx 
import copy
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism
import itertools

def add_labeled_edge(labeled_graph, node1, node2):
    cluster_attrs = nx.get_node_attributes(labeled_graph, 'cluster')
    flag_attrs = nx.get_node_attributes(labeled_graph, 'flag')
    node_attrs_list = [(cluster_attrs[node1], flag_attrs[node1]), (cluster_attrs[node2], flag_attrs[node2])]
    labeled_graph.add_edge(node1, node2, node_attrs=set(node_attrs_list))

def equiv_nodes(node1, node2):
    return ((node1['flag'], node1['cluster']) == (node2['flag'], node2['cluster']))

def equiv_edges(edge1, edge2):
    return (edge1['node_attrs'] == edge2['node_attrs'])

def label_is_isomorphic(graph1, graph2):
    return nx.is_isomorphic(graph1, graph2, node_match=equiv_nodes,
        edge_match=equiv_edges)

def compose_graphs_on_flag(graph1, graph2, flag_iso):
    # Check that the given graphs have distinct node names
    assert len(set(list(graph1.nodes())).intersection(set(list(graph2.nodes())))) == 0
    renamed_graph2 = nx.relabel_nodes(graph2, flag_iso)
    return nx.compose(graph1, renamed_graph2)

def combination(n, k):
    product = 1
    for i in range(1, k+1):
        product *= n - i + 1
        product //= i
    return product

# In the form of a list of pairs (SingleGraph, float coefficient)
def draw_graph_form(graph_form):
    if len(graph_form) < 4:
        num_cols = len(graph_form)
        num_rows = 1
    else:
        num_cols = int(math.sqrt(len(graph_form)))
        num_rows = int(math.ceil(float(len(graph_form)) / num_cols))
    f, axs = plt.subplots(num_rows, num_cols, sharey=True)
    
    for idx in range(len(graph_form)):
        single_graph, coeff = graph_form[idx]
        if len(graph_form) == 1:
            ax = axs
        elif len(graph_form) < 4:
            ax = axs[idx]
        else:
            ax = axs[idx // num_cols, idx % num_cols]
        single_graph.draw(extra_text=('%.5g' % coeff), ax=ax, draw=False)
    plt.subplots_adjust(hspace=0.7)
    plt.show()
    plt.close()

class SingleGraph:
    def __init__(self, orig_graph, cluster_dict, flag_dict):
        # Empty the graph attributes from orig_graph
        self.graph = nx.Graph()
        self.graph.add_nodes_from(list(orig_graph.nodes()))
        self.graph.add_edges_from(list(orig_graph.edges()))

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
            add_labeled_edge(self.graph, n1, n2)

    def draw(self, flag_size=400, node_size=200, extra_text=None, ax=None, draw=True):
        cluster_labels = nx.get_node_attributes(self.graph, 'cluster')
        cluster_labels = { k: (cluster_labels[k] + 1) for k in cluster_labels }
        color_map = { 1: 'r', 2: 'gold', 3: 'cornflowerblue' }
        cluster_colors = [color_map[cluster_labels[k]] for k in sorted(cluster_labels.keys())]

        pos = nx.circular_layout(self.graph)
        nx.draw_networkx_nodes(self.graph, pos, node_size=node_size, 
            node_color=cluster_colors, ax=ax)

        flag_status = nx.get_node_attributes(self.graph, 'flag')
        flags = [n for n in flag_status if flag_status[n]]
        nx.draw_networkx_nodes(self.graph, pos, nodelist=flags, 
            node_shape='d', node_color=cluster_colors, node_size=flag_size, ax=ax)
        nx.draw_networkx_labels(self.graph, pos, labels=cluster_labels, ax=ax)

        nx.draw_networkx_edges(self.graph, pos, ax=ax)

        if extra_text is not None:
            ax.set_title(extra_text)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
        if draw:
            plt.show()
            plt.close()
        return ax

    def is_isomorphic(self, network2):
        # Both nodes and edges need to have equivalence properties
        return label_is_isomorphic(self.graph, network2.graph)

    # Outputs a renamed version of graph2
    def rename_disjoint_graphs(self, single_graph2):
        max_node1 = max(list(self.graph.nodes())) - min(list(single_graph2.graph.nodes()))
        # Rename the nonflag graph2 nodes so they're all distinct from those in graph1
        renaming = { n: (n + max_node1 + 1) for n in list(single_graph2.graph.nodes()) }
        renamed_single_graph2 = copy.deepcopy(single_graph2)
        renamed_single_graph2.graph = nx.relabel_nodes(single_graph2.graph, renaming)
        return renamed_single_graph2

    def get_flag(self):
        flag_attributes = nx.get_node_attributes(self.graph, 'flag')
        node_list = [node for node in flag_attributes if flag_attributes[node]]
        return self.graph.subgraph(node_list)

    def multiply(self, single_graph2_orig):
        # Ensure none of the nodes have the same name
        single_graph2 = self.rename_disjoint_graphs(single_graph2_orig)
        # Check the flags are isomorphic
        flag1 = self.get_flag()
        flag2 = single_graph2.get_flag()
        # important I put flag2 first, so that the isomorphism maps from flag2 to flag1
        flag_matcher = isomorphism.GraphMatcher(flag2, flag1)  
        if not flag_matcher.is_isomorphic():
            print('Flag with nodes {} and flag with nodes {} not isomorphic.'.format(
                list(flag1.nodes()), list(flag2.nodes())))
            raise AssertionError

        non_flag_nodes1 = list(set(self.graph.nodes()) - set(flag1.nodes()))
        non_flag_nodes2 = list(set(single_graph2.graph.nodes()) - set(flag2.nodes()))

        # Compile all possible graphs of size 
        possible_combined_graphs = []

        # Check among all isomorphisms between the two flags 
        for flag_iso in flag_matcher.isomorphisms_iter():
            empty_combined_graph = compose_graphs_on_flag(
                self.graph, single_graph2.graph, flag_iso)

            # Add any set of edges that you want also, to make a graph with nonzero coefficient
            for edge_bools in itertools.product(range(2), repeat=len(non_flag_nodes1) * len(non_flag_nodes2)):
                combined_graph = copy.deepcopy(empty_combined_graph)
                cnt = 0
                for node1 in non_flag_nodes1:
                    for node2 in non_flag_nodes2:
                        if edge_bools[cnt]:
                            # Don't forget to add the edge attributes!
                            add_labeled_edge(combined_graph, node1, node2)
                        cnt += 1

                # but still double check that there's no isomorphisms we've already seen
                is_included = False
                for prev_graph in possible_combined_graphs:
                    if label_is_isomorphic(prev_graph, combined_graph):
                        is_included = True
                        break
                if not is_included:
                    possible_combined_graphs.append(combined_graph)

        graph_product = []

        # Find the probabilities/coefficients for each graph
        #total_combo_count = combination(len(non_flag_nodes1 + non_flag_nodes2), len(non_flag_nodes1))
        for combined_graph in possible_combined_graphs:
            iso_count = 0
            total_combo_count = 0
            for cand1_nodes_tup in itertools.combinations(non_flag_nodes1 + non_flag_nodes2, len(non_flag_nodes1)):
                cand1_nodes = list(cand1_nodes_tup)
                cand2_nodes = list(set(non_flag_nodes1 + non_flag_nodes2) - set(cand1_nodes))
                cand1 = nx.subgraph(combined_graph, cand1_nodes + list(flag1.nodes()))
                cand2 = nx.subgraph(combined_graph, cand2_nodes + list(flag1.nodes()))
                
                # Check if the clusters are the same in cand1 & graph1, same with cand2 & graph2
                cand1_clusters = sorted(nx.get_node_attributes(cand1, 'cluster').values()) 
                orig1_clusters = sorted(nx.get_node_attributes(self.graph, 'cluster').values()) 
                cand2_clusters = sorted(nx.get_node_attributes(cand2, 'cluster').values()) 
                orig2_clusters = sorted(nx.get_node_attributes(single_graph2.graph, 'cluster').values()) 
                if cand1_clusters == orig1_clusters and cand2_clusters == orig2_clusters:
                    total_combo_count += 1

                if label_is_isomorphic(cand1, self.graph) and label_is_isomorphic(cand2, single_graph2.graph):
                    iso_count += 1
            cluster_attrs = nx.get_node_attributes(combined_graph, 'cluster')
            flag_attrs = nx.get_node_attributes(combined_graph, 'flag')
            combined_singlegraph = SingleGraph(combined_graph, cluster_attrs, flag_attrs)
            graph_product.append((combined_singlegraph, float(iso_count) / total_combo_count))
        return graph_product

    def check_cluster_independence(self):
        cluster_labels = nx.get_node_attributes(self.graph, 'cluster')
        for n1, n2 in self.graph.edges():
            if cluster_labels[n1] == cluster_labels[n2]:
                return False 
        return True

    # TODO: test this!
    def count_labeled_edges(self, label1, label2):
        count = 0
        edge_attrs = nx.get_edge_attributes(self.graph, 'node_attrs')
        for edge in edge_attrs:
            if set([attr[0] for attr in edge_attrs[edge]]) == set([label1, label2]):
                count += 1
        return count

    def count_labeled_pairs(self, label1, label2):
        cluster_labels = nx.get_node_attributes(self.graph, 'cluster')
        if label1 == label2:
            cnt = len([node for node in cluster_labels if cluster_labels[node] == label1])
            return cnt * (cnt - 1) // 2
        else:
            cnt1 = len([node for node in cluster_labels if cluster_labels[node] == label1])
            cnt2 = len([node for node in cluster_labels if cluster_labels[node] == label2])
            return cnt1 * cnt2

    def unflag_graph(self):
        cluster_labels = nx.get_node_attributes(self.graph, 'cluster')
        flag_labels = { idx: False for idx in self.graph.nodes() }
        return SingleGraph(self.graph, cluster_labels, flag_labels)

    # NOTE: this is a hack to averaging down a flag of size 1
    def average_single_flag(self):
        # Assert that there is only one flag node
        flag_attrs = nx.get_node_attributes(self.graph, 'flag')
        orig_flag = [idx for idx in flag_attrs if flag_attrs[idx]][0]
        if sum(flag_attrs.values()) != 1:
            raise AssertionError('Flag is not of size 1')

        unflagged = self.unflag_graph()
        # Iterate thru all nodes
        iso_count = 0
        total_count = 0
        cluster_labels = nx.get_node_attributes(unflagged.graph, 'cluster')
            

        for flag in unflagged.graph.nodes():
            if cluster_labels[flag] != cluster_labels[orig_flag]:
                continue
            total_count += 1
            # For each one, make a new graph with that node as the flag, then test isomorphism
            flag_labels = { idx: False for idx in unflagged.graph.nodes() }
            flag_labels[flag] = True
            new_graph = SingleGraph(unflagged.graph, cluster_labels, flag_labels)
            if self.is_isomorphic(new_graph):
                iso_count += 1

        return (unflagged, float(iso_count) / total_count)


    #def average_flags(self):
        # probability is only among labels that match correctly (so don't account for when sigma is mapped to wrong cluster labels)
        # unlike in multiplication, we care about distinct isomorphisms of a subset (ordering of vertices)

# TODO: check if all the above functions work with a nonsymmetric flag (i.e. w/ different isomorphisms)
# TODO: incorporate into a testing structure rather than a single function
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

    print('\n*******\n')


    # Test composing/multiplying graphs with nonisomorphic flags
    print 'Muliplying/composing graphs with wrong flags: ',
    try:
        point_raw_graph = nx.Graph()
        point_raw_graph.add_node(1)
        point_graph = SingleGraph(point_raw_graph, {1: 1}, { 1: False })
        point_graph.multiply(v_graph)
        print('no exception thrown (bad!)')
    except AssertionError:
        print('exception thrown (good!)')

    # TODO: test the renaming + composing functions

    # test the flag_notes examples here (both with and without flags)
    double_raw1 = nx.Graph()
    double_raw1.add_nodes_from([1, 2])
    double_raw2 = copy.deepcopy(double_raw1)
    double_raw2.add_edge(1, 2)
    double_clusters = { 1: 1, 2: 1 }
    double_no_flags = { 1: False, 2: False }
    double_no_flag_graph1 = SingleGraph(double_raw1, double_clusters, double_no_flags)
    double_no_flag_graph2 = SingleGraph(double_raw2, double_clusters, double_no_flags)
    double_no_flag_product = double_no_flag_graph1.multiply(double_no_flag_graph2)
    draw_graph_form(double_no_flag_product)

    double_flags = { 1: True, 2: False }
    double_flag_graph1 = SingleGraph(double_raw1, double_clusters, double_flags)
    double_flag_graph2 = SingleGraph(double_raw2, double_clusters, double_flags)
    double_flag_product_01 = double_flag_graph1.multiply(double_flag_graph2)
    draw_graph_form(double_flag_product_01)

    double_flag_product_11 = double_flag_graph2.multiply(double_flag_graph2)
    draw_graph_form(double_flag_product_11)

    # TODO: add tests for flags larger than 1 vertex
    # TODO: add tests to ensure constructor is emptying the attributes

    draw_graph_form([v_graph.average_single_flag()]) # should have coefficient 1

    uniform_v_raw_graph = nx.Graph()
    uniform_v_raw_graph.add_nodes_from([1, 2, 3])
    uniform_v_raw_graph.add_edges_from([(1, 2), (1, 3)])
    uniform_cluster_dict = { 1: 1, 2: 1, 3: 1 }
    uniform_flag_dict = { 1: True, 2: False, 3: False }
    uniform_v_graph = SingleGraph(uniform_v_raw_graph, uniform_cluster_dict, uniform_flag_dict)
    draw_graph_form([uniform_v_graph.average_single_flag()])

    uniform_flag_dict = { 1: False, 2: True, 3: False }
    uniform_v_graph = SingleGraph(uniform_v_raw_graph, uniform_cluster_dict, uniform_flag_dict)
    draw_graph_form([uniform_v_graph.average_single_flag()])

if __name__ == '__main__':
    run_tests()

