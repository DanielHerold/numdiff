import numpy as np
import matplotlib.pyplot as plt

values = np.loadtxt('netzwerk3.data')
time = values[:,0]
alpha = values[:,1]
alpha_dot = values[:,2]

plt.plot(time, alpha)
plt.show()

plt.plot(time, alpha_dot)
plt.show()
