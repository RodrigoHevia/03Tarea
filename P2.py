'''
El siguiente script resuelve y plotea las ecuaciones de Lorenz para el caso del
Atractor de Lorenz utilizando el metodo de Runge-Kutta de orden 4
'''


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode

fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

def Lorenz(x, y, z, sigma = 10, rho = 28, beta = 8./3) :

    dx = sigma*(y - x)
    dy = x*(rho - z) - y
    dz = x*y-beta*z

    return dx, dy, dz

x0 = 0
y0 = 0
z0 = 0
f0 = [x0, y0, z0]

sol = ode(Lorenz, f0).set_integrator('dopri5').set_initial_value(f0)

print "hola mundo", sol

#ax.plot(sol)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
