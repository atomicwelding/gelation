# Bethe lattice

Bethe lattice has been a convenient tool to describe gelation. It is a special case of general bond-percolation, namely Bernouilli percolation model. We model the Bethe lattice as a graph, with the following building rules, 

* Each monomer is represented by a node.
* We define functionality of the monomer to be the maximum number of links a node can possess, and note it $f$.
* The root node is considered to be linked to a parent node which is not represented, so it can have a maximum of $f-1$ children.
* Each link has a probability $p \in [0,1]$ to occur.
* A monomer can only bond with the $f$ neighbouring monomers.
* Graph is assumed to be infinite.

From that, analytical and numerical results can be derived. We define $p_c$ the critical probability above which ($p > p_c$) a large network spans.


## Analytical results 

In construction. 


## Simulation [file : `bethe.py`]

We represent the lattice as an adjacency matrix of dimension $n_\text{monomers}$ and construct it following building rules described above.


### Average degree of polymerization

This quantity is pretty straightforward to obtain : we count the number of nodes in each connected components and divide by the number of components.

### Pgel, Psol

This quantity is the ratio between the degree of polymerization of the gel and the total number of monomers. $P_\text{sol} = 1 - P_\text{gel}$.

### Conductance

We propose a model to compute the conductance of the system. Once percolation occurs, we focus on the largest connected component of the network, referred to as the *gel*. Specifically, we randomly select a leaf node (a node with degree 1) from the gel, which we designate as the *root* node. Next, we introduce a new node, referred to as the *superleaf*, into the gel graph. This superleaf is connected to all other leaf nodes in the gel except for the root node. To compute the conductance, we apply a voltage difference of 1V between the root and superleaf nodes, solving Kirchhoff’s laws for the resulting network. In practice, this involves calculating the resistance distance using the Moore-Penrose pseudo-inverse formalism.

This model aims to approximate the experimental situation by representing the gel as a conductive network, where current flows through the largest connected cluster. The addition of the superleaf is a way to simulate boundary conditions, similar to how electrodes might interact with multiple points in a real system. On the other hand, the choice of the root node is made random because the exact location of current injection in real systems is often not well-defined or uniform. Additionally, the topology of the resistance network—specifically, which resistances are in parallel or in series—can vary significantly depending on the selected root node. By making this choice random and relying on statistical averaging across multiple simulations, we aim to reduce potential biases and ensure that the computed conductance reflects the overall network behavior rather than being overly influenced by specific configurations.

### Results [file : `building_statistics.py`, `plotting.py`]

We perform statistics for different $f$ with $n_\text{monomers} = 5000$ monomers. For each system, we're taking 200 $p$ from $0$ to $1$ evenly distributed. We're averaging over 100 such systems.
