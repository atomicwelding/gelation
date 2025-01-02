import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

# parameters
nnodes = 3000 
lattice = np.zeros((nnodes, nnodes)) 

p = 1 # probability of bonding
f = 4  # functionnaliy of each monomer


shouldConnect = lambda: np.random.choice([0, 1], p=[1 - p, p])

# build the lattice 
for n in range(nnodes):
    for i in range(1, f):  # each monomer (node) connects to f - 1 children
        child = (f - 1) * n + i # child index 
        if child < nnodes:  # ensure child index is within bounds
            connect = shouldConnect()
            lattice[n, child] = connect
            lattice[child, n] = connect  

# visualize the graph
G = nx.from_numpy_array(lattice)
plt.figure(figsize=(12, 8))
nx.draw(G, node_size=10, with_labels=False)
plt.show()
