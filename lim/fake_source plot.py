import matplotlib.pyplot as plt
import matplotlib
import numpy as np

x = np.linspace(0,53.3, 127)

y = np.zeros(len(x))
for i in range(0, len(x-1)):
    if x[i] >= 0 and x[i] <= 9.68:
        y[i] = 5e12
    elif x[i] > 9.68 and x[i] <= 17.99:
        y[i] = 3.5e12
    elif x[i] >= 25.02 and x[i] <= 37.06:
         y[i] = 5e12
    else:
        y[i] = 0

matplotlib.rc('font', size= 25)
plt.scatter(x, y, color = 'b')
plt.title("Interface Profile Suggested by CPs", fontsize = 27, fontweight = 'bold')
plt.xlabel("Distance from target [m]", fontsize = 27, fontweight = 'bold')
plt.ylabel("Number of Impurity Ions", fontsize = 27, fontweight = 'bold')

plt.show()

