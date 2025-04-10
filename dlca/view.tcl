# Load the XYZ trajectory
mol new data.xyz type xyz waitfor all

# Supprime la représentation par défaut
mol delrep 0 top

# Ajoute une représentation VDW colorée par nom
mol representation VDW 0.3 12.0
mol color Name
mol selection all
mol material Opaque
mol addrep top

# Remet la vue par défaut et centre sur la boîte
display resetview
display projection Orthographic
scale to 1.0

# Fond noir pour bon contraste
color Display Background black

# Définir manuellement la boîte périodique (10x10x10)
pbc set {10 10 10} -all
pbc box
pbc box -shiftcenter {+0.5 +0.5 +0.5}
pbc display box