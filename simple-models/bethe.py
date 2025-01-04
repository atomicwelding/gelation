import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


# parameters
nmonomers = 1000
lattice = np.zeros((nmonomers, nmonomers)) 

p = 1. # probability of bonding
f = 3 # functionnaliy of each monomer


shouldConnect = lambda: np.random.choice([0, 1], p=[1 - p, p])

# build the lattice
print("=== BUILDING LATTICE ===")
for n in range(nmonomers):
    for i in range(1, f):  # each monomer (node) connects to f - 1 children
        child = (f - 1) * n + i # child index 
        if child < nmonomers:  # ensure child index is within bounds
            connect = shouldConnect()
            lattice[n, child] = connect
            lattice[child, n] = connect
    if(n % 1000 == 0):
        print(n/nmonomers * 100, "% step")


print("=== MAKING STATS OF THE GRAPH ===") 
# create a graph
G = nx.from_numpy_array(lattice)

# statistics
components = nx.connected_components(G)
polymers = [component for component in components]
degrees_of_polymerization = [len(polymer) for polymer in polymers]

avg_degree = np.mean(degrees_of_polymerization)
max_degree = max(degrees_of_polymerization)

# gel and conductance calculation
gel_index = degrees_of_polymerization.index(max_degree)
gel = G.subgraph(polymers[gel_index]).copy()


# we choose randomly a leaf as the root
leaves = [node for node in gel.nodes() if gel.degree(node) == 1]
root = np.random.choice(leaves)
leaves.remove(root)

# we connect a "superleaf" to every other leaves 
super_leaf = "Superleaf"
gel.add_node(super_leaf)

for leaf in leaves:
    gel.add_edge(leaf, super_leaf)


# we compute the resistivity
print("test")
resistivity = nx.resistance_distance(gel, root, super_leaf)


# fractions
Pgel = max_degree / nmonomers
Psol = 1 - Pgel

# print("=== RESULTS ===")
print(f"{nmonomers = }\t\t{p = :.2f}")
print(f"{Psol = :.2f}\t\t{Pgel = :.2f}")
print(f"avg deg = {avg_degree:.2f}\t\tmax deg = {max_degree}")
print(f"conductance = { 1/resistivity:.2f}")
if len(sys.argv) > 1 and sys.argv[1] == "-v":
    # visualize the graph
    plt.figure(figsize=(12, 8))
    nx.draw(G, node_size=10, with_labels=False)
    plt.show()
