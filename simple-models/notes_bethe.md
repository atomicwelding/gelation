# Bethe lattice

Bethe lattice has been a convenient tool to describe gelation. It is a special case of general bond-percolation, namely Bernouilli percolation model. We model the Bethe lattice as a graph, with the following building rules, 

* Each monomer is represented by a node.
* We define functionality of the monomer to be the maximum number of links a node can possess, and note it $f$.
* The root node is considered to be linked to a parent node which is not represented, so it can have a maximum of $f-1$ children.
* Each link has a probability $p \in [0,1] $ to occur.
* A monomer can only bond with the $f$ neighbouring monomers.
* Graph is assumed to be infinite.

From that, analytical and numerical results can be derived. We define $p_c$ the critical probability above which ($p > p_c$) a large network spans.


## Analytical results 

In construction. 


## Numerical results [file : `bethe.py`]

We represent the lattice as an adjacency matrix of dimension $n_\text{nodes}$ and construct it following building rules described above.