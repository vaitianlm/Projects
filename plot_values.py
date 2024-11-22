import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 16})


N = 28
T = np.zeros(N)
eps_mean = np.zeros(N)
m_abs_mean = np.zeros(N)
heat_cap = np.zeros(N)
mag_sus = np.zeros(N)
eps_sq = np.zeros(N)
m_sq = np.zeros(N)

L = np.array([40, 60, 80, 100])
color = ['orange', 'blue', 'green', 'purple']

for j in range(len(L)):
    with open(f"L{L[j]}_N200000_28_values.txt") as infile:
        i=0
        for line in infile:
            T[i] = line.split()[0]
            eps_mean[i] = line.split()[1]
            m_abs_mean[i] = line.split()[2]
            heat_cap[i] = line.split()[3]
            mag_sus[i] = line.split()[4]
            eps_sq[i] = line.split()[5]
            m_sq[i] = line.split()[6]
            i+=1

    plt.plot(T, mag_sus, 'o', label = f"$L={L[j]}$", color = color[j])
    
    fit = np.polyfit(T, mag_sus, 2)
    fitfunc = fit[0]*T**2+ fit[1]*T + fit [2]
    top = fitfunc[np.argmax(fitfunc)]
    T_top = T[np.argmax(fitfunc)]
    print(T_top)
    plt.plot(T, fit[0]*T**2+ fit[1]*T + fit [2], linestyle = "dashed", color = color[j])
    plt.plot(T_top, top, 'x', color = "red", markersize = 14)


plt.legend()
plt.grid()
plt.xlabel("$T$ [$\,J/k_B$]")
plt.ylabel("$\\lange |m| \rangle$ [$\,J^{-1}$]")
plt.subplots_adjust(left = 0.167, bottom = 0.138)
plt.legend(fontsize = 16)
plt.savefig("fig_magsus_critical.pdf")
plt.show()
