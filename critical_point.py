import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import linregress
plt.rcParams.update({"font.size": 16})



T_c = np.array([2.28, 2.29, 2.30, 2.32])
T_c = [2.2907, 2.2967, 2.2981, 2.3141]
L = np.array([100, 80, 60, 40])

linreg = linregress(1/L, T_c)

print(linreg.intercept)


plt.plot(1/L, T_c, 'o')

plt.grid()
L = np.array([999990, 100, 80, 60, 40])
plt.plot(1/L, linreg.slope/L + linreg.intercept, linestyle = 'dashed')
plt.plot(0, linreg.intercept, 'x', color = "black", label = "Extrapolation", markersize = 10)
plt.plot(0, 2.269, 'x', label = "Analytical", color = "r", markersize = 10)
plt.xlabel('$1/L$')
plt.ylabel('$T_C(L^{-1})$ $[\,J/k_B]$')
plt.subplots_adjust(left = 0.160, bottom = 0.138)
plt.legend()
plt.savefig("fig_critical.pdf")
plt.show()
