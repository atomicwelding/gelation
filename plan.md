# Rapport — Percolation, Transition Sol-Gel et Modélisation DLCA

Une ébauche pour mon plan de rapport. j'ai aussi quelques notes dans mon cahier de manip.

## Introduction

Les transitions sol-gel sont au cœur de nombreux systèmes physiques (colloïdes, polymères, verres). La percolation permet d’interpréter cette transition comme l’émergence d’une structure macroscopique connectée. Deux approches sont complémentaires : des modèles théoriques (réseau de Bethe) et des simulations (DLCA). Ce rapport met en lumière ces approches et leur lien.

---

## 1. Approche théorique : comprendre la percolation

### 1.1 Exemple simple : percolation sur réseau carré
- Introduction à la percolation site/bond
- Définition du seuil critique, notion de cluster infini

### 1.2 Réseau de Bethe
- Avantages analytiques (absence de boucle)
- Calcul du seuil de percolation
- Analyse de la connectivité et exponent critique β

### 1.3 Conductivité et analogie avec un réseau de résistances
- Transition de phase → transition de conductivité
- Loi de puissance autour du seuil : σ ∝ (p - p_c)^t

---

## 2. Approche numérique : simulation DLCA

### 2.1 DLCA — Diffusion-Limited Cluster Aggregation
- Description du modèle (agrégation + diffusion)
- Passage d’un réseau statique à une dynamique réelle

### 2.2 Du réseau de Bethe à DLCA
- Bethe donne la connectivité
- DLCA y ajoute : espace, temps, collisions
- Une manière de **physiquer** un modèle abstrait

### 2.3 Article de Gimel et al.
- Présentation de leurs hypothèses et résultats
- Problèmes méthodologiques identifiés :
  - Absence de critère de convergence clair
  - Temps d’arrêt arbitraire ?

### 2.4 Implémentation Fortran
- Description générale (pas à pas)
- Stratégies de convergence (test de croissance nulle)
- Optimisations introduites (critère d'arrêt dynamique, flush I/O)

### 2.5 Ouverture vers la dynamique moléculaire
- DLCA comme amorce d’un gel réaliste
- Possibilité de transférer vers GROMACS / LAMMPS

---

## 3. Résultats

### 3.1 P(φ₀, L) — Mesure de la transition
- Graphe : probabilité de gel vs densité initiale
- Comparaison entre tailles L
- Mise en évidence d’un seuil φ_c
- Comparaison avec les résultats de Gimel

### 3.2 Temps de gel
- Mesure de t_gel via 1/Nc
- Distribution des temps de gel
- Réflexion sur la dépendance en L et φ₀
- Évocation de données expérimentales si disponibles

### 3.3 Analyse complémentaire
- Histogramme des tailles de clusters
- Estimation de la dimension fractale ?
- Nombre de clusters vs temps

---

## Discussion

- Limites du modèle : bruit, convergence lente
- Différences notables avec l’article de Gimel
- La question du **temps d’arrêt** comme paramètre critique
- Importance de définir un critère reproductible

---

## Conclusion

- DLCA complète les modèles de percolation statique
- La simulation révèle des subtilités du gelage non visibles dans Bethe
- Perspectives : vers la modélisation réaliste de gels (MD)
