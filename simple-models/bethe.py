import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

# parameters
nmonomers = 10000
lattice = np.zeros((nmonomers, nmonomers)) 

p = 0.55 # probability of bonding
f = 3 # functionnaliy of each monomer


shouldConnect = lambda: np.random.choice([0, 1], p=[1 - p, p])

# build the lattice 
for n in range(nmonomers):
    for i in range(1, f):  # each monomer (node) connects to f - 1 children
        child = (f - 1) * n + i # child index 
        if child < nmonomers:  # ensure child index is within bounds
            connect = shouldConnect()
            lattice[n, child] = connect
            lattice[child, n] = connect  

# create a graph
G = nx.from_numpy_array(lattice)

# statistics
polymers = list(nx.connected_components(G))
degrees_of_polymerization = [len(polymer) for polymer in polymers]

avg_degree = np.mean(degrees_of_polymerization)
max_degree = max(degrees_of_polymerization)

Pgel = max_degree / nmonomers
Psol = 1 - Pgel

print(f"{Psol = :.2f}\t{Pgel = :.2f}")
print(f"avg deg = {avg_degree:.2f}")
if len(sys.argv) > 1 and sys.argv[1] == "-v":
    # visualize the graph
    plt.figure(figsize=(12, 8))
    nx.draw(G, node_size=10, with_labels=False)
    plt.show()
