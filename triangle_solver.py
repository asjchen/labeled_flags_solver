# Semidefinite program for triangle case

import argparse
import networkx as nx
from numpy import matrix
import Irene
from graph_utils import SingleGraph
import itertools

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

#def flag_alg_constraints(vector_indices):
    # take any flag, then compute all the possible pairs of flag vectors
    # multiply all those pairs, apply averaging operator, correspond the coefficients with the vector entries

# given an edge, the sum of the graphs containing that edge
def rho_density_constraints(vector_indices, rho):
    target_edges = [(0, 1), (1, 2), (0, 2)]
    lhs_constraints = []
    rhs_constraints = []
    for label1, label2 in target_edges:
        valid_indices = set()
        for idx in range(len(vector_indices)):
            if vector_indices[idx].contains_labeled_edge(label1, label2):
                valid_indices.add(idx)
        lhs_constraints.append([matrix([[0]]) for j in range(len(vector_indices))])
        for idx in valid_indices:
            lhs_constraints[-1][idx] = matrix([[-1]])
        rhs_constraints.append(matrix([[-1 * rho]]))
    return lhs_constraints, rhs_constraints


def main():
    parser = argparse.ArgumentParser()
    # TODO: add description to this
    parser.add_argument('rho', type=float)
    args = parser.parse_args()

    # compute the indexing
    vector_indices = generate_indexing()
    print len(vector_indices) # should be 56

    # obtain constraints via functions, each of which should return a tuple, which will get separated
    nonneg_constraints = nonnegative_entry_constraints(vector_indices)
    print 'Nonnegative Constraints: {}'.format(len(nonneg_constraints[0])) # should be 56

    label_sum_constraints = label_sums_to_one(vector_indices)
    print 'Label Sum Constraints (summing to 1): {}'.format(len(label_sum_constraints[0])) # should be 20 (10 for each inequality direction)

    cluster_ind_constraints = cluster_independence_constraints(vector_indices)
    print 'Cluster Independendence Constraints: {}'.format(len(cluster_ind_constraints[0])) # should be 27

    rho_constraints = rho_density_constraints(vector_indices, args.rho)
    print 'Rho Density Constraints: {}'.format(len(rho_constraints[0])) # should be 3

    # flag algebra constraints
    # convert constraints to Irene form
    # write the objective
    # solve dat!



if __name__ == '__main__':
    main()
