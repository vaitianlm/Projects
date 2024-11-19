import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 16})


N = 10000
n_MC_cycles = np.zeros(N)
eps_T10 = np.zeros(N)
mean_eps_T10 = np.zeros(N)

eps_T24 = np.zeros(N)
mean_eps_T24 = np.zeros(N)


with open("T2.4_L20_u_samples.txt") as infile:
    i=0
    for line in infile:
        n_MC_cycles[i] = line.split()[0]
        eps_T10[i] = line.split()[1]
        mean_eps_T10[i] = line.split()[2]
        i+=1

# with open("distribution_T24_N1000000.txt") as infile:
#     i=0
#     for line in infile:
#         eps_T24[i] = line.split()[1]
#         mean_eps_T24[i] = line.split()[2]
#         i+=1


plt.plot(n_MC_cycles, eps_T10, '-', color='#377eb8', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, mean_eps_T10, '-', linewidth=2.0, color='#377eb8', label='$T=1.0$ $J/k_{B}$')

# plt.plot(n_MC_cycles, eps_T24, '-', color='#e41a1c', alpha=0.4, linewidth=1.0)
# plt.plot(n_MC_cycles, mean_eps_T24, '-', linewidth=2.0, color='#e41a1c', label='$T=2.4$ $J/k_{B}$')

plt.subplots_adjust(left = 0.167, bottom = 0.138)
plt.xlabel("MCMC-cycles done")
plt.ylabel("Lattice energy per spin [$J$]") # Ask about dimention of epsilon
plt.grid()
plt.legend(fontsize = 13.5)
plt.xscale("log")
# plt.savefig("epsilon_T1.pdf")
plt.show()

plt.style.use("seaborn")
bins = np.arange(np.min(eps_T24), np.max(eps_T24), 0.009999999999)
print(bins[:20])

plt.hist(eps_T24, density = True, histtype = "stepfilled", bins = bins)
plt.savefig("epsilon_disr_T1.pdf")
plt.show()

bins = np.arange(np.min(eps_T10), np.max(eps_T10), 0.00099999999)
print(bins[:20])

plt.hist(eps_T10, density = True, histtype = "stepfilled", bins = bins)
# plt.savefig("epsilon_disr_T24.pdf")
plt.show()