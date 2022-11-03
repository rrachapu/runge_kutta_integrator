import numpy as np
import rk_ode
import matplotlib.pyplot as plt

print("hello")
def spring(t, x, props):
    k = props[0]
    c = props[1]
    x_eq = props[2]
    
    x_pos = x[0]
    x_vel = x[1]
    
    dv = -1*k*(x_pos - x_eq) -1*c*x_vel
    dx = x_vel + dv

    return np.array([dx,dv])
    
t = 0
x0 = np.array([0,.2])
t = np.linspace(0,10,500)
dt = t[1] - t[0]
x = np.zeros((t.shape[0], 2))
i = 1
eps = 1e-3
props = [2, .1, 0]
while (i < t.shape[0]):
    print(i)
    x[i] = rk_ode.ode45(2, x[i-1], t[i-1], t[i], eps, dt, 1e-3, spring, props)
    i += 1

plt.figure(0)
plt.plot(t, np.transpose(x[:,0]))

