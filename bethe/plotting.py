import numpy as np
import matplotlib.pyplot as plt

fp_conductance = "statistics_conductance.txt"
fp_avg = "statistics_avg_polymerization.txt"
fp_Pgel = "statistics_Pgel.txt"

N_TEST, N_POINTS, nmonomers, f = (0,0,0,0)
with open(fp_conductance) as file:
    N_TEST, N_POINTS, nmonomers, f = map(int, file.readline().strip().split())


# plotting conductance G(p)
data_conductance = np.loadtxt(fp_conductance, skiprows=1)
avg_conductance  = np.mean(data_conductance[:,1:], axis=1)

plt.plot(data_conductance[250:,0], avg_conductance[250:], 'k')
plt.xlabel(r"$p$")
plt.ylabel(r"$G(p)$")
plt.show()

# plotting avg polymerization
data_avg = np.loadtxt(fp_avg, skiprows=1)
avg_avg = np.mean(data_avg[:,1:], axis=1)

plt.plot(data_avg[:,0], avg_avg, 'k')
plt.xlabel(r"$p$")
plt.ylabel("avg deg poly")
plt.show()

# plotting Pgel
data_Pgel = np.loadtxt(fp_Pgel, skiprows=1)
avg_Pgel = np.mean(data_Pgel[:,1:], axis=1)

plt.plot(data_Pgel[:,0], avg_Pgel, 'k')
plt.xlabel(r"$p$")
plt.ylabel(r"$P_\text{gel}$")
plt.show()



