import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

N_TEST = 100
N_POINTS = 200


# parameters
nmonomers = 5000
f = 3

lattice = np.zeros((nmonomers, nmonomers))


conductances = np.zeros((N_TEST, N_POINTS))
avgs_polymerization = np.zeros((N_TEST, N_POINTS))
Pgels = np.zeros((N_TEST,N_POINTS))

ps = np.linspace(0,1,num=N_POINTS)



# main loop
print("=== simulating ... ===")
shouldConnect = lambda: np.random.choice([0, 1], p=[1 - p, p])
for test in range(N_TEST):
    for idx, p in enumerate(ps):
        
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
        if not leaves:
            conductances[test,idx] = 0
            continue
        
        root = np.random.choice(leaves)
        leaves.remove(root)

        # we connect a "superleaf" to every other leaves 
        super_leaf = "Superleaf"
        gel.add_node(super_leaf)
        
        for leaf in leaves:
            gel.add_edge(leaf, super_leaf)


        # we compute the resistivity
        resistivity = nx.resistance_distance(gel, root, super_leaf)


        # fractions
        Pgel = max_degree / nmonomers

        Pgels[test,idx] = Pgel
        avgs_polymerization[test,idx] = avg_degree
        conductances[test, idx] = 1/resistivity


print("=== writing data ...===")
with open("statistics_conductance.txt","w") as saving:
    saving.write(f"{N_TEST} {N_POINTS} {nmonomers} {f}\n")
    for idx, p in enumerate(ps):
        saving.write(f"{p} ")
        for test in range(N_TEST):
            saving.write(f"{conductances[test,idx]} ")
        saving.write("\n")

with open("statistics_avg_polymerization.txt", "w") as saving:
    saving.write(f"{N_TEST} {N_POINTS} {nmonomers} {f}\n")
    for idx, p in enumerate(ps):
        saving.write(f"{p} ")
        for test in range(N_TEST):
            saving.write(f"{avgs_polymerization[test,idx]} ")
        saving.write("\n")
            
with open("statistics_Pgel.txt", "w") as saving:
    saving.write(f"{N_TEST} {N_POINTS} {nmonomers} {f}\n")
    for idx, p in enumerate(ps):
        saving.write(f"{p} ")
        for test in range(N_TEST):
            saving.write(f"{Pgels[test,idx]} ")
        saving.write("\n")           
