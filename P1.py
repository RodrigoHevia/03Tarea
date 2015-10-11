#! /usr/bin/env python

'''
El siguiente script resuelve la ecuacion diferencial del oscilador de Van der Pol
luego de un cambio de variable, mediante el metodo de Runge Kutta de orden 3.
'''

import matplotlib.pyplot as plt
import numpy as np

plt.figure(1)
plt.clf()

def rk3_2(f, s0, sf, y0, v0, n) :

    '''
    Esta funcion resuelve la ecuacion diferencial de segundo orden "f" mediante
    el metodo de Runge-Kutta de orden 3, con condiciones iniciales "y0" y "v0"
    en el intervalo [s0,sf] con resolucion "n"
    '''

    s = np.linspace(s0, sf, n, endpoint = True)
    y = np.zeros(n)
    v = np.zeros(n)

    h = (sf - s0)/n

    y[0] = y0
    v[0] = v0

    for i in range(1, n) :

        k1 = h*v[i-1]
        l1 = h*f(s[i-1], y[i-1], v[i-1])
        k2 = h*(v[i-1] + 0.5*l1)
        l2 = h*f(s[i-1] + 0.5*h, y[i-1] + 0.5*k1, v[i-1] + 0.5*l1)
        k3 = h*(v[i-1] + 0.5*l2)

        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3)/5
        v[i] = v[i-1] + (l1 + 2*l2)/3

    return s, y

def f(s, y, v) :

    u2 = 1.699
    dv = - y - u2 * ( y**2 - 1) * v

    return dv

s, y = rk3_2(f, 0, 20*np.pi, 0.1, 0, 1000) # dice t no s

plt.plot(s, y, label='Oscilador de Van der Pol')

plt.xlabel('s')
plt.ylabel('y')
plt.legend()

plt.draw()
plt.show()
