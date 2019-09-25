# Semidefinite program for triangle case

import argparse
import networkx as nx
from numpy import matrix, zeros, sqrt
from numpy.linalg import norm, eigvals
import Irene
from graph_utils import SingleGraph, draw_graph_form
import itertools
import matplotlib.pyplot as plt

def combination(n, k):
    prod = 1
    for i in range(1, k + 1):
        prod *= n + 1 - i
        prod /= i
    return prod

# edge_bools[i] denotes whether there's an edge from i to i+1
def make_triangle_graph(clusters, edge_bools, flags=[]):
    raw_graph = nx.Graph()
    raw_graph.add_nodes_from(range(3))
    raw_graph.add_edges_from([(i, (i + 1) % 3) for i in range(3) if edge_bools[i]])
    cluster_dict = { i: clusters[i] for i in range(3) }
    flag_dict = { i: False for i in range(3) }
    flag_dict.update({ i: True for i in flags })
    return SingleGraph(raw_graph, cluster_dict, flag_dict)

# construct all possible graphs and generate an indexing
def generate_indexing():
    # Iterate over all possible singlegraphs
    vector_order = []
    for clusters in itertools.product(range(3), repeat=3):
        for edge_bools in itertools.product(range(2), repeat=3):
            triangle_graph = make_triangle_graph(clusters, edge_bools)
            is_new = True
            for prev_graph in vector_order:
                if prev_graph.is_isomorphic(triangle_graph):
                    is_new = False
                    break
            if is_new:
                vector_order.append(triangle_graph)
    return vector_order
    # TODO: can we do fast lookups indexed by graph? like some hash combined with a dictionary?
    # Or fast based on number of edges + cluster names, then linear search?

# target should be a SingleGraph object
def vector_entry_lookup(vector_indices, target):
    # search thru all until we see one that's isomorphic
    for i in range(len(vector_indices)):
        if target.is_isomorphic(vector_indices[i]):
            return i
    return None # TODO: potentially throw an error here

def visualize_vector(vector_indices, num_vector, tolerance=1e-8):
    # check every index, and if its entry is nonzero, add to list
    graph_form = []
    for idx in range(len(vector_indices)):
        if num_vector[idx] > tolerance:
            graph_form.append((vector_indices[idx], num_vector[idx]))
    draw_graph_form(graph_form)

def visualize_indexing(vector_indices):
    visualize_vector(vector_indices, [1] * len(vector_indices))

def visualize_constraint(vector_indices, lhs_constraint):
    visualize_vector(vector_indices, [norm(m) for m in lhs_constraint])

def nonnegative_entry_constraints(vector_indices):
    lhs_constraints = []
    rhs_constraints = []
    for i in range(len(vector_indices)):
        # put one-hot vectors on LHS and 0 on RHS
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        lhs_constraints[i][i] = matrix([[1]])
        rhs_constraints.append(matrix([[0]]))
    return lhs_constraints, rhs_constraints

def label_sums_to_one(vector_indices):
    # since we're working with triangles, let's enumerate cluster label configs manually
    cluster_configs = [(i, j, j) for j in range(3) for i in range(3)]
    cluster_configs += [(0, 1, 2)]
    lhs_constraints = []
    rhs_constraints = []
    for clusters in cluster_configs:
        valid_indices = set()
        # iterate thru all possible edge configs, see what its index is
        for edge_bools in itertools.product(range(2), repeat=3):
            # on LHS: all zeros, ones where the labels are correct
            # on RHS: one
            triangle_graph = make_triangle_graph(clusters, edge_bools)
            idx = vector_entry_lookup(vector_indices, triangle_graph)
            if idx is None:
                raise AssertionError(('Could not find graph with clusters {} '
                    'and edge_bools {} in vector'.format(clusters, edge_bools)))
            valid_indices.add(idx)
        # Inequality in the first direction
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        for idx in valid_indices:
            lhs_constraints[-1][idx] = matrix([[1]])
        rhs_constraints.append(matrix([[1]]))
        # Inequality in the other direction (to ensure equality constraint)
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        for idx in valid_indices:
            lhs_constraints[-1][idx] = matrix([[-1]])
        rhs_constraints.append(matrix([[-1]]))
    return lhs_constraints, rhs_constraints

# every graph with an edge within an edge should have zero count
def cluster_independence_constraints(vector_indices):
    lhs_constraints = []
    rhs_constraints = []
    for i in range(len(vector_indices)):
        if not vector_indices[i].check_cluster_independence():
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            lhs_constraints[-1][i] = matrix([[-1]])
            rhs_constraints.append(matrix([[0]]))
    return lhs_constraints, rhs_constraints

def make_pair_graph(flag_cluster, partner_cluster, edge_bool):
    raw = nx.Graph()
    raw.add_nodes_from([1, 2])
    if edge_bool:
        raw.add_edge(1, 2)
    clusters = { 1: flag_cluster, 2: partner_cluster }
    flags = { 1: True, 2: False }
    return SingleGraph(raw, clusters, flags)

# Note: in K3 case, this doesn't accomplish any further constraints, as there 
# exists a soln on the remaining constraints that also happens to satisfy this
# set of flag constraints
def flag_alg_constraints(vector_indices):
    lhs_constraints = []
    rhs_constraints = []
    partners = [(cluster, edge) for edge in range(2) for cluster in range(3)]
    for flag_cluster in range(3):
        lhs_constraints.append([matrix(zeros((6, 6))) for _ in range(len(vector_indices))])
        rhs_constraints.append(matrix(zeros((6, 6))))
        for idx1 in range(len(partners)):
            cluster1, edge1 = partners[idx1]
            pair1 = make_pair_graph(flag_cluster, cluster1, edge1)
            for idx2 in range(len(partners)):
                cluster2, edge2 = partners[idx2]
                pair2 = make_pair_graph(flag_cluster, cluster2, edge2)
                product = pair1.multiply(pair2)
                # multiply the two coefficients 
                for sing_graph, mult_coeff in product:
                    avg_graph, avg_coeff = sing_graph.average_single_flag()
                    graph_idx = vector_entry_lookup(vector_indices, avg_graph)
                    lhs_constraints[flag_cluster][graph_idx][idx1, idx2] += mult_coeff * avg_coeff
                    # print lhs_constraints[flag_cluster][graph_idx]
    return lhs_constraints, rhs_constraints


# given an edge, the sum of the graphs containing that edge
def rho_density_constraints_old(vector_indices, rho):
    target_edges = [(0, 1), (1, 2), (0, 2)]
    cluster_configs = [(i, j, j) for j in range(3) for i in range(3)]
    cluster_configs += [(0, 1, 2)]
    lhs_constraints = []
    rhs_constraints = []

    for label1, label2 in target_edges:
        for clusters in cluster_configs:
            empty_graph = make_triangle_graph(clusters, [0, 0, 0])
            num_target_pairs = empty_graph.count_labeled_pairs(label1, label2)
            rho_dict = {}
            for edge_bools in itertools.product(range(2), repeat=3):
                triangle_graph = make_triangle_graph(clusters, edge_bools)
                num_target_edges = triangle_graph.count_labeled_edges(label1, label2)
                if num_target_edges not in rho_dict:
                    rho_dict[num_target_edges] = set()
                graph_idx = vector_entry_lookup(vector_indices, triangle_graph)
                rho_dict[num_target_edges].add(graph_idx)
            print label1, label2
            print clusters
            for num_target_edges in rho_dict:
                thresh = (rho ** num_target_edges) * ((1 - rho) ** (num_target_pairs - num_target_edges))
                thresh *= combination(num_target_pairs, num_target_edges)
                
                print num_target_edges, thresh

                lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
                for idx in rho_dict[num_target_edges]:
                    lhs_constraints[-1][idx] = matrix([[1]])
                rhs_constraints.append(matrix([[thresh]]))

                lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
                for idx in rho_dict[num_target_edges]:
                    lhs_constraints[-1][idx] = matrix([[-1]])
                rhs_constraints.append(matrix([[-thresh]]))

        # rho_dict = {}
        
        # for clusters in cluster_configs:
        #     for edge_bools in itertools.product(range(2), repeat=3):
        #         triangle_graph = make_triangle_graph(clusters, edge_bools)
        #         num_target_pairs = triangle_graph.count_labeled_pairs(label1, label2)
        #         num_target_edges = triangle_graph.count_labeled_edges(label1, label2)
        #         num_target_nonedges = num_target_pairs - num_target_edges
        #         if (num_target_edges, num_target_nonedges) not in rho_dict:
        #             rho_dict[(num_target_edges, num_target_nonedges)] = []
        #         graph_idx = vector_entry_lookup(vector_indices, triangle_graph)
        #         rho_dict[(num_target_edges, num_target_nonedges)].append(graph_idx)

        # for num_target_edges, num_target_nonedges in rho_dict:
        #     # Enforce equality to rho^p * (1-rho)^q
        #     lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        #     for idx in rho_dict[(num_target_edges, num_target_nonedges)]:
        #         lhs_constraints[-1][idx] = matrix([[1]])
        #     rhs_constraints.append(matrix([[(rho ** num_target_edges) * ((1 - rho) ** num_target_nonedges)]]))

        #     lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        #     for idx in rho_dict[(num_target_edges, num_target_nonedges)]:
        #         lhs_constraints[-1][idx] = -matrix([[1]])
        #     rhs_constraints.append(matrix([[-(rho ** num_target_edges) * ((1 - rho) ** num_target_nonedges)]]))
        #     # make the constraints

    # for each pair of clusters
    #   for each cluster configuration
    #       count the number of possible edges k
    #       for edge_bools among the target edges:
    #           compute the number of existing edges m
    #           compute kCm * (1-p)^(k-m) * p^m as RHS (exact?)
    #           for remaining_edge_bools among the rest of the edges:
    #               add to the constraint


    # for label1, label2 in target_edges:
    #     for clusters in cluster_configs:
    #         valid_indices = set()
    #         for edge_bools in itertools.product(range(2), repeat=3):
            
    #             triangle_graph = make_triangle_graph(clusters, edge_bools)
    #             idx = vector_entry_lookup(vector_indices, triangle_graph)
    #             if idx is None:
    #                 raise AssertionError(('Could not find graph with clusters {} '
    #                     'and edge_bools {} in vector'.format(clusters, edge_bools)))    
    #             if vector_indices[idx].contains_labeled_edge(label1, label2):
    #                 valid_indices.add(idx)
    #         lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
    #         for idx in valid_indices:
    #             lhs_constraints[-1][idx] = matrix([[-1]])
    #         rhs_constraints.append(matrix([[-1 * rho]]))
    return lhs_constraints, rhs_constraints


# NOTE: currently very specifically for the triangle case...
def rho_density_constraints(vector_indices, rho):
    target_edges = [(0, 1), (1, 2), (0, 2)]
    lhs_constraints = []
    rhs_constraints = []

    # Ignore the 0-0-0 case (all three vertices from same cluster)

    # Examine 0-0-1 case
    for orig_label1, orig_label2 in target_edges:
        for label1, label2 in [(orig_label1, orig_label2), (orig_label2, orig_label1)]:
            clusters = (label1, label2, label2)
            # 2 edges
            two_graphs = [make_triangle_graph(clusters, [1, i, 1]) for i in range(2)]
            two_indices = [vector_entry_lookup(vector_indices, two_graph) for two_graph in two_graphs]
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for two_idx in two_indices:
                lhs_constraints[-1][two_idx] = matrix([[1]])
            rhs_constraints.append(matrix([[rho * rho]]))
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for two_idx in two_indices:
                lhs_constraints[-1][two_idx] = matrix([[-1]])
            rhs_constraints.append(matrix([[-rho]]))

            # 1 edge
            one_graphs = [make_triangle_graph(clusters, [1, i, 0]) for i in range(2)]
            one_indices = [vector_entry_lookup(vector_indices, one_graph) for one_graph in one_graphs]
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for one_idx in one_indices:
                lhs_constraints[-1][one_idx] = matrix([[1]])
            rhs_constraints.append(matrix([[0]]))
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for one_idx in one_indices:
                lhs_constraints[-1][one_idx] = matrix([[-1]])
            rhs_constraints.append(matrix([[-2 * rho * (1 - rho)]]))

            # 0 edges
            zero_graphs = [make_triangle_graph(clusters, [0, i, 0]) for i in range(2)]
            zero_indices = [vector_entry_lookup(vector_indices, zero_graph) for zero_graph in zero_graphs]
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for zero_idx in zero_indices:
                lhs_constraints[-1][zero_idx] = matrix([[1]])
            rhs_constraints.append(matrix([[(1 - rho) ** 2]]))
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for zero_idx in zero_indices:
                lhs_constraints[-1][zero_idx] = matrix([[-1]])
            rhs_constraints.append(matrix([[-(1 - rho)]]))

            # extra constraint?
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for zero_idx in zero_indices:
                lhs_constraints[-1][zero_idx] = matrix([[1]])
            for one_idx in one_indices:
                lhs_constraints[-1][one_idx] = matrix([[2]])
            for two_idx in two_indices:
                lhs_constraints[-1][two_idx] = matrix([[1]])
            rhs_constraints.append(matrix([[1]]))
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for zero_idx in zero_indices:
                lhs_constraints[-1][zero_idx] = matrix([[-1]])
            for one_idx in one_indices:
                lhs_constraints[-1][one_idx] = matrix([[-2]])
            for two_idx in two_indices:
                lhs_constraints[-1][two_idx] = matrix([[-1]])
            rhs_constraints.append(matrix([[-1]]))

            # another constraint?
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for zero_idx in zero_indices:
                lhs_constraints[-1][zero_idx] = matrix([[1]])
            for one_idx in one_indices:
                lhs_constraints[-1][one_idx] = matrix([[1]])
            rhs_constraints.append(matrix([[1 - rho]]))
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for zero_idx in zero_indices:
                lhs_constraints[-1][zero_idx] = matrix([[-1]])
            for one_idx in one_indices:
                lhs_constraints[-1][one_idx] = matrix([[-1]])
            rhs_constraints.append(matrix([[-1 + rho]]))

    # Finally the 0-1-2 case
    # Use fact we have at most one edge existing in this configuration
    perms = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
    for clusters in perms:
        # if the first edge is the target edge

        # Empty edge
        zero_indices = []
        for edge_bools in itertools.product(range(2), repeat=2):
            zero_graph = make_triangle_graph(clusters, (0,) + edge_bools)
            zero_indices.append(vector_entry_lookup(vector_indices, zero_graph))
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        for zero_idx in zero_indices:
            lhs_constraints[-1][zero_idx] = matrix([[1]])
        rhs_constraints.append(matrix([[1 - rho]]))
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        for zero_idx in zero_indices:
            lhs_constraints[-1][zero_idx] = matrix([[-1]])
        rhs_constraints.append(matrix([[-(1 - rho)]]))

        # Nonempty edge
        one_indices = []
        for edge_bools in itertools.product(range(2), repeat=2):
            one_graph = make_triangle_graph(clusters, (1,) + edge_bools)
            one_indices.append(vector_entry_lookup(vector_indices, one_graph))
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        for one_idx in one_indices:
            lhs_constraints[-1][one_idx] = matrix([[1]])
        rhs_constraints.append(matrix([[rho]]))
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        for one_idx in one_indices:
            lhs_constraints[-1][one_idx] = matrix([[-1]])
        rhs_constraints.append(matrix([[-rho]]))

    return lhs_constraints, rhs_constraints


def solve_triangle_problem(rho, verbose=False):
    # Compute the indexing
    vector_indices = generate_indexing()
    print len(vector_indices) # should be 56

    # obtain constraints via functions, each of which should return a tuple, which will get separated
    nonneg_constraints = nonnegative_entry_constraints(vector_indices)
    if verbose:
        print 'Nonnegative Constraints: {}'.format(len(nonneg_constraints[0])) # should be 56

    label_sum_constraints = label_sums_to_one(vector_indices)
    if verbose:
        print 'Label Sum Constraints (summing to 1): {}'.format(len(label_sum_constraints[0])) # should be 20 (10 for each inequality direction)

    cluster_ind_constraints = cluster_independence_constraints(vector_indices)
    if verbose:
        print 'Cluster Independendence Constraints: {}'.format(len(cluster_ind_constraints[0])) # should be 27

    flag_constraints = flag_alg_constraints(vector_indices)
    if verbose:
        print 'Flag Constraints: {}'.format(len(flag_constraints[0])) # should be 3

    rho_constraints = rho_density_constraints(vector_indices, rho)
    if verbose:
       print 'Rho Density Constraints: {}'.format(len(rho_constraints[0])) # should be 30 (3 for each cluster configuration)

    lhs_constraints = nonneg_constraints[0] + label_sum_constraints[0] + cluster_ind_constraints[0] + rho_constraints[0] + flag_constraints[0]
    rhs_constraints = nonneg_constraints[1] + label_sum_constraints[1] + cluster_ind_constraints[1] + rho_constraints[1] + flag_constraints[1]


    # convert constraints to Irene form
    constraint_blocks = []
    for i in range(len(vector_indices)):
        constraint_blocks.append([])
        for j in range(len(lhs_constraints)):
            constraint_blocks[i].append(lhs_constraints[j][i])
    SDP = Irene.sdp('cvxopt')
    SDP.Option('feastol', 1e-9)
    SDP.Option('maxiters', 1000)
    SDP.AddConstantBlock(rhs_constraints)
    for i in range(len(vector_indices)):
        SDP.AddConstraintBlock(constraint_blocks[i])

    # write the objective: minimizing independent sets
    obj_raw_graph = nx.Graph()
    obj_raw_graph.add_nodes_from(range(3))
    obj_clusters = { i: i for i in range(3) }
    obj_flags = { i: False for i in range(3) }
    obj_graph = SingleGraph(obj_raw_graph, obj_clusters, obj_flags)
    obj_idx = vector_entry_lookup(vector_indices, obj_graph)
    obj = [0] * len(vector_indices)
    obj[obj_idx] = 1
    SDP.SetObjective(obj)


    # y = [ 1.00000000e+00,  4.06614762e-12,  3.38194173e-12, -6.52881171e-13,
    #     6.66256585e-01,  5.66092511e-04,  3.33177323e-01, -5.83826685e-12,
    #    -5.56071968e-12, -5.77429640e-12,  6.66256585e-01,  5.66092511e-04,
    #     3.33177323e-01, -5.83824722e-12, -5.56070062e-12, -5.77427691e-12,
    #     6.66256585e-01,  5.66092511e-04, -5.83827484e-12, -5.56072745e-12,
    #     3.33177323e-01, -5.77430434e-12, -2.41884097e-11,  3.33224259e-01,
    #     3.33224259e-01,  9.04762200e-05,  3.33224259e-01,  9.04762200e-05,
    #     9.04762200e-05,  5.57953073e-05,  6.66256585e-01,  5.66092511e-04,
    #    -5.83827941e-12, -5.56073189e-12,  3.33177323e-01, -5.77430888e-12,
    #     1.00000000e+00,  4.06615984e-12,  3.38196175e-12, -6.52859719e-13,
    #     6.66256585e-01,  5.66092511e-04,  3.33177323e-01, -5.83827519e-12,
    #    -5.56072778e-12, -5.77430469e-12,  6.66256585e-01,  5.66092511e-04,
    #    -5.83826627e-12, -5.56071912e-12,  3.33177323e-01, -5.77429583e-12,
    #     1.00000000e+00,  4.06615837e-12,  3.38195933e-12, -6.52862302e-13]
    # flag1 = sum([y[i] * flag_constraints[0][0][i] for i in range(56)])
    # print flag1
    # print '********'
    # print eigvals(flag1)
    # print '********'


    # solve dat!
    SDP.solve()
    if verbose:
        print SDP.Info
        visualize_vector(vector_indices, SDP.Info['y'])
    if 'PObj' in SDP.Info:
        return SDP.Info['PObj']
    return None


def main():
    # parser = argparse.ArgumentParser()
    # # TODO: add description to this
    # parser.add_argument('rho', type=float)
    # args = parser.parse_args()

    # pair_raw_graph = nx.Graph()
    # pair_raw_graph.add_nodes_from([1, 2])
    # pair1 = SingleGraph(pair_raw_graph, {1: 0, 2: 0}, {1: True, 2: False})
    # pair2 = SingleGraph(pair_raw_graph, {1: 0, 2: 1}, {1: True, 2: False})
    # product = pair1.multiply(pair2)
    # draw_graph_form(product)
    # second = product[1][0]
    # draw_graph_form([second.average_single_flag()])
    


    intervals = 100
    orig_rhos = [i / float(intervals) for i in range(intervals + 1)]

    #rhos = [2 - (1 + sqrt(5.0)) / 2]
    orig_rhos = [0.33333333]

    rhos = []
    pobjs = []
    for rho in orig_rhos:
        #pobj = solve_triangle_problem(rho, verbose=False)
        pobj = solve_triangle_problem(rho, verbose=True)
        if pobj is not None:
            rhos.append(rho)
            pobjs.append(pobj)
        print rho
        print pobjs[-1]
    plt.plot(rhos, pobjs)
    plt.show()
    plt.close()



    


if __name__ == '__main__':
    main()
