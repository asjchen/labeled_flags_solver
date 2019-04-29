# Labeled Flag Algebra Solver for Multipartite Turan Problem

## Python Requirements
* cvxopt
* Irene
* networkx

## Todo List:
* Flag Algebras
    * decide consistent and nonredundant representation of graphs and flags (might be good to use NetworkX (or whatever I used in 224W) or another package)
    * find a way to visually display the graphs and flags
    * write a test suite
    * multiply two flag graphs
    * average a flag graph
* Semidefinite program
    * decide indexing of the decision vector
    * write constraints for nonnegative entries
    * write constraints for summing similar labels to 1
    * write constraints for semidefinite matrices
    * write constraints for edge densities
    * write objective
    * solve the problem!
