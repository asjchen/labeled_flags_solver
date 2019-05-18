# Semidefinite program for triangle case

import argparse
import networkx as nx
from numpy import matrix, zeros
from numpy.linalg import norm
import Irene
from graph_utils import SingleGraph, draw_graph_form
import itertools
import matplotlib.pyplot as plt

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
    return lhs_constraints, rhs_constraints


# TODO: this is written incorrectly -- what exactly needs to sum to rho?

# given an edge, the sum of the graphs containing that edge
def rho_density_constraints(vector_indices, rho):
    target_edges = [(0, 1), (1, 2), (0, 2)]
    cluster_configs = [(i, j, j) for j in range(3) for i in range(3)]
    cluster_configs += [(0, 1, 2)]
    lhs_constraints = []
    rhs_constraints = []

    for label1, label2 in target_edges:
        for clusters in cluster_configs:
            valid_indices = set()
            for edge_bools in itertools.product(range(2), repeat=3):
            
                triangle_graph = make_triangle_graph(clusters, edge_bools)
                idx = vector_entry_lookup(vector_indices, triangle_graph)
                if idx is None:
                    raise AssertionError(('Could not find graph with clusters {} '
                        'and edge_bools {} in vector'.format(clusters, edge_bools)))    
                if vector_indices[idx].contains_labeled_edge(label1, label2):
                    valid_indices.add(idx)
            lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
            for idx in valid_indices:
                lhs_constraints[-1][idx] = matrix([[-1]])
            rhs_constraints.append(matrix([[-1 * rho]]))
    return lhs_constraints, rhs_constraints

def solve_triangle_problem(rho, verbose=False):
    # Compute the indexing
    vector_indices = generate_indexing()
    print len(vector_indices) # should be 56

    #visualize_indexing(vector_indices)

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

    # solve dat!
    SDP.solve()
    if verbose:
        print SDP.Info
    visualize_vector(vector_indices, SDP.Info['y'])
    return SDP.Info['PObj']


def main():
    # parser = argparse.ArgumentParser()
    # # TODO: add description to this
    # parser.add_argument('rho', type=float)
    # args = parser.parse_args()

    intervals = 1000
    rhos = [i / float(intervals) for i in range(intervals + 1)]

    rhos = [0.334]

    pobjs = []
    for rho in rhos:
        pobjs.append(solve_triangle_problem(rho, verbose=True))
        print rho
        print pobjs[-1]
    plt.plot(rhos, pobjs)
    plt.show()
    plt.close()



    


if __name__ == '__main__':
    main()
