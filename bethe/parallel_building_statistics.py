import sys
import networkx as nx
import numpy as np
import time
from multiprocessing import Pool, cpu_count

MAX_CORE = cpu_count()


# === PARAMÈTRES GÉNÉRAUX ===
N_TEST = 100      # Nombre de tests pour chaque probabilité
N_POINTS = 100       # Nombre de points de probabilité
nmonomers = 10000     # Nombre total de monomères
f = 6               # Fonctionnalité de chaque monomère
ps = np.linspace(0, 1, num=N_POINTS)  # Liste des probabilités


# === FONCTIONS UTILITAIRES ===
def build_lattice(nmonomers, f, p):
    """Construit un graphe aléatoire basé sur une probabilité de liaison p."""
    lattice = np.zeros((nmonomers, nmonomers))
    should_connect = lambda: np.random.choice([0, 1], p=[1 - p, p])
    
    for n in range(nmonomers):
        for i in range(1, f):  # Chaque monomère (nœud) connecte à f-1 enfants
            child = (f - 1) * n + i
            if child < nmonomers:  # Vérifie si l'index est dans les limites
                connect = should_connect()
                lattice[n, child] = connect
                lattice[child, n] = connect
    return lattice


def analyze_gel(G):
    """Identifie le gel et calcule ses propriétés."""
    components = nx.connected_components(G)
    polymers = [component for component in components]
    degrees_of_polymerization = [len(polymer) for polymer in polymers]
    
    max_degree = max(degrees_of_polymerization)
    gel_index = degrees_of_polymerization.index(max_degree)
    gel = G.subgraph(polymers[gel_index]).copy()
    
    return gel, max_degree, np.mean(degrees_of_polymerization)


def calculate_conductance(gel, nmonomers):
    """Calcule la conductance du gel."""
    leaves = [node for node in gel.nodes if gel.degree[node] == 1]
    if not leaves:
        return 0, 0  # Aucun leaf, conductance impossible à calculer
    
    root = np.random.choice(leaves)
    leaves.remove(root)
    super_leaf = "Superleaf"
    gel.add_node(super_leaf)
    
    for leaf in leaves:
        gel.add_edge(leaf, super_leaf)
    
    try:
        resistivity = nx.resistance_distance(gel, root, super_leaf)
        conductance = 1 / resistivity
    except ZeroDivisionError:
        conductance = 0
    
    Pgel = len(gel.nodes) / nmonomers
    return conductance, Pgel


def simulate_for_p(p):
    print("new p ...")
    """Effectue la simulation pour une probabilité donnée p."""
    conductances = []
    avgs_polymerization = []
    Pgels = []

    for _ in range(N_TEST):
        # Construire le réseau
        lattice = build_lattice(nmonomers, f, p)
        G = nx.from_numpy_array(lattice)
        
        # Analyser le gel
        gel, max_degree, avg_degree = analyze_gel(G)
        
        # Calculer la conductance et Pgel
        conductance, Pgel = calculate_conductance(gel, nmonomers)
        
        # Stocker les résultats
        conductances.append(conductance)
        avgs_polymerization.append(avg_degree)
        Pgels.append(Pgel)

    return p, conductances, avgs_polymerization, Pgels


# === MAIN ===
if __name__ == '__main__':
    start_simulation = time.time()
    print("=== Simulating in parallel ===")

    # Exécution parallèle
    with Pool(processes=MAX_CORE) as pool:
        results = pool.map(simulate_for_p, ps)

    # Récupération des résultats
    conductances = np.zeros((N_TEST, N_POINTS))
    avgs_polymerization = np.zeros((N_TEST, N_POINTS))
    Pgels = np.zeros((N_TEST, N_POINTS))

    for idx, (p, cond, avg_poly, pgel) in enumerate(results):
        conductances[:, idx] = cond
        avgs_polymerization[:, idx] = avg_poly
        Pgels[:, idx] = pgel

    end_simulation = time.time()
    print(f"Simulation completed in {end_simulation - start_simulation:.2f}s")

    # === ÉCRITURE DES RÉSULTATS ===
    start_writing = time.time()
    print("=== Writing data ===")

    with open("statistics_conductance.txt", "w") as saving:
        saving.write(f"{N_TEST} {N_POINTS} {nmonomers} {f}\n")
        for idx, p in enumerate(ps):
            saving.write(f"{p} " + " ".join(map(str, conductances[:, idx])) + "\n")

    with open("statistics_avg_polymerization.txt", "w") as saving:
        saving.write(f"{N_TEST} {N_POINTS} {nmonomers} {f}\n")
        for idx, p in enumerate(ps):
            saving.write(f"{p} " + " ".join(map(str, avgs_polymerization[:, idx])) + "\n")

    with open("statistics_Pgel.txt", "w") as saving:
        saving.write(f"{N_TEST} {N_POINTS} {nmonomers} {f}\n")
        for idx, p in enumerate(ps):
            saving.write(f"{p} " + " ".join(map(str, Pgels[:, idx])) + "\n")

    end_writing = time.time()
    print(f"Writing completed in {end_writing - start_writing:.2f}s")
    print(f"Total time: {end_writing - start_simulation:.2f}s")
