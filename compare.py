import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 16})

N = 100000

n_MC_cycles = np.zeros(N)
eps = np.zeros(N)
eps_mean = np.zeros(N)
m_abs_mean = np.zeros(N)
heat_cap = np.zeros(N)
mag_sus = np.zeros(N)



with open(f"T1_L2_N100000_o_samples.txt") as infile:
    i=0
    for line in infile:
        n_MC_cycles[i] = line.split()[0]
        eps[i] = line.split()[1]
        eps_mean[i] = line.split()[2]
        m_abs_mean[i] = line.split()[3]
        heat_cap[i] = line.split()[4]
        mag_sus[i] = line.split()[5]
        i+=1

z = 4*np.cosh(8)+12
eps_anal = -2*np.sinh(8)/(np.cosh(8)+3)
m_abs_anal = (2*np.exp(8)+4)/z

C_V_anal = (64 * (3 * z - 16)) / (z**2)
mag_sus_anal = 16* (3+3*np.exp(8)+np.exp(-8))/z**2

print('relative error |m|%: ', (abs((m_abs_anal-m_abs_mean[20000])/m_abs_anal))*100)
print('Correct |m|: ', m_abs_anal)
print('Wrong |m|: ', m_abs_mean[20000])


plt.plot(n_MC_cycles[:20000], mag_sus[:20000], label = "Metropolis algorithm")
plt.axhline(mag_sus_anal, color='red', linestyle='dashed', label = "Analytical")
plt.legend()
plt.grid()
plt.xlabel("MCMC-cycles done")
plt.ylabel("$C_V(T=1)/N_{spins}$ [$k_B$]")
plt.subplots_adjust(left = 0.205, bottom = 0.138)
plt.legend(fontsize = 16)
# plt.savefig("fig_heatcap_compare.pdf")
# plt.xscale("log")
plt.show()