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

def Lorenz_Atractor(t, f, sigma = 10., rho = 28., beta = 8./3) :

    '''
    Sistema de ecuaciones diferenciales del Atractor de Lorenz
    '''

    dx = sigma*(f[1] - f[0])
    dy = f[0]*(rho - f[2]) - f[1]
    dz = f[0]*f[1]-beta*f[2]

    return [dx, dy, dz]

def Jac(t, f, sigma = 10., rho = 28., beta = 8./3) :

    '''
    Jacobiano para el sistema de ecuaciones del Atractor de Lorenz
    '''

    jac = [[-sigma, sifma, 0], [rho-f[2], -1, -f[0]], [f[1], f[0], -beta]]

    return jac

x0 = 31.5
y0 = 25.4
z0 = 12.1

f0 = [x0, y0, z0]
t0 = 0

sol = ode(Lorenz_Atractor, Jac).set_integrator('dopri5')
sol.set_initial_value(f0, t0)

dt = 0.001
tf = 30.0

x = np.zeros((tf-t0)/dt)
y = np.zeros((tf-t0)/dt)
z = np.zeros((tf-t0)/dt)

i = 0
while sol.successful() and sol.t < tf :

    x[i], y[i], z[i] = sol.integrate(sol.t+dt)
    i += 1

ax.plot(x, y, z, label = "Atractor de Lorenz" )

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.title('Atractor de Lorenz')
plt.show()

fig.savefig('Atractor de Lorenz.jpg')
