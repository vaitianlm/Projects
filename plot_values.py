import numpy as np 
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 16})


N = 31
T = np.zeros(N)
eps_mean = np.zeros(N)
m_abs_mean = np.zeros(N)
heat_cap = np.zeros(N)
mag_sus = np.zeros(N)
eps_sq = np.zeros(N)
m_sq = np.zeros(N)



with open("L80_N10000000_31_values.txt") as infile:
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

plt.plot(T, mag_sus, '-')
plt.show()