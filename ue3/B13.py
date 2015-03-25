import numpy as np
import matplotlib.pyplot as plt

values = np.loadtxt('B13.out')

time=values[:,0]
y = values[:,1]

plt.plot(time, y)
plt.show()

