import numpy as np
import matplotlib.pyplot as plt

values = np.loadtxt('pendel_mp.data')
time = values[:,0]
alpha = values[:,1]
alpha_dot = values[:,2]

plt.plot(time, alpha)
plt.show()

q=alpha
p=alpha_dot
H=np.cos(q)*(-1)+(p**2)/(2)

plt.plot(time,H)
plt.show()

plt.plot(q,p)
plt.show()

