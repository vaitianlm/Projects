import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 16})


N = 10000
n_MC_cycles = np.zeros(N)
eps_T10_o = np.zeros(N)
mean_eps_T10_o = np.zeros(N)

eps_T10_u = np.zeros(N)
mean_eps_T10_u = np.zeros(N)

eps_T24_o = np.zeros(N)
mean_eps_T24_o = np.zeros(N)

eps_T24_u = np.zeros(N)
mean_eps_T24_u = np.zeros(N)

with open("burn_in_T1_ordered.txt") as infile:
    i=0
    for line in infile:
        n_MC_cycles[i] = line.split()[0]
        eps_T10_o[i] = line.split()[1]
        mean_eps_T10_o[i] = line.split()[2]
        i+=1

with open("burn_in_T1_unordered.txt") as infile:
    i=0
    for line in infile:
        eps_T10_u[i] = line.split()[1]
        mean_eps_T10_u[i] = line.split()[2]
        i+=1

with open("burn_in_T2.4_ordered.txt") as infile:
    i=0
    for line in infile:
        eps_T24_o[i] = line.split()[1]
        mean_eps_T24_o[i] = line.split()[2]
        i+=1

with open("burn_in_T2.4_unordered.txt") as infile:
    i=0
    for line in infile:
        eps_T24_u[i] = line.split()[1]
        mean_eps_T24_u[i] = line.split()[2]
        i+=1

plt.plot(n_MC_cycles, eps_T10_u, '-', color='#377eb8', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, eps_T10_o, '-', color='#4daf4a', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, eps_T24_u, '-', color='#e41a1c', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, eps_T24_o, '-', color='#984ea3', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, mean_eps_T10_u, '-', linewidth=2.0, color='#377eb8', label='$T=1.0$ $J/k_{B}$, unordered')
plt.plot(n_MC_cycles, mean_eps_T10_o, '-', linewidth=2.0, color='#4daf4a', label='$T=1.0$ $J/k_{B}$, ordered')
plt.plot(n_MC_cycles, mean_eps_T24_u, '-', linewidth=2.0, color='#e41a1c', label='$T=2.4$ $J/k_{B}$, unordered')
plt.plot(n_MC_cycles, mean_eps_T24_o, '-', linewidth=2.0, color='#984ea3', label='$T=2.4$ $J/k_{B}$, ordered')

plt.subplots_adjust(left = 0.167, bottom = 0.138)
plt.xlabel("MCMC-cycles done")
plt.ylabel("Lattice energy per spin [$J$]") # Ask about dimention of epsilon
plt.grid()
plt.legend(fontsize = 13.5)
plt.xscale("log")

plt.savefig("fig_burn_in.pdf")
plt.show()