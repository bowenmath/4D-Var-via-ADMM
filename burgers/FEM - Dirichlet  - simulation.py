#------------------------------------------------------
# This program uses finitie element method to simulate
# Burgers equation.
#------------------------------------------------------
# By Bowen Li
# October 7, 2024
#------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
import time


def g(x, t):

   """brugers model"""
   # Setting up vector
   d = np.zeros(N)
   # Loops over indices (with operations and Python underflow indexing handling edge cases)
   d[0] = 0
   d[N - 1] = 0
   for i in range(1, N-1):
      d[i] = (1/6) * (x[i-1]**2 + x[i]*x[i-1] - x[i+1]*x[i] - x[i+1]**2) + (nu/(dx)) * (x[i-1] - 2*x[i] + x[i+1])
   return np.matmul(G, d)

# The iteration map x_{i} = f(x_{i-1})
def f(x, dt):
   t = 0
   return x + dt*g(x, t)


time_start = time.time()

# These are our constants
N = 101  # Number of variables

nu = 0.005

dt = 0.01
num_steps = 200
dx = np.pi/(N-1)

G = np.zeros([N, N])
for i in range(1, N-1):
   G[i, i-1] = dx/6
   G[i, i] = 2*dx/3
   G[i, i+1] = dx/6
G[0, 0] = 2*dx/3
G[0, 1] = dx/6
G[N-1, N-2] = dx/6
G[N-1, N-1] = 2*dx/3

G = np.linalg.inv(G)

x = np.linspace(0, np.pi, N)
x0 = np.sin(x)  # Initial state (equilibrium)
#x0 += np.random.normal(0, 1, N)  # Add small perturbation to the first variable

'''t_span = [0.0, num_steps * dt]
#t = np.arange(0.0, (num_steps + 1)*dt, dt)
t = np.linspace(0.0, num_steps * dt, num_steps + 1)
x = np.linspace(-np.pi, np.pi, num_steps + 1)

xs = solve_ivp(g, t_span, x0, t_eval=t).y
xs = xs.T'''

# Time integration
#x = np.linspace(0, np.pi, N)
xs = np.zeros([num_steps+1, N])
xs[0] = x0
for i in range(num_steps):
   xs[i+1] = f(xs[i], dt)

print(xs[0])


plt.figure(0)
plt.plot(x, xs[0], label='Initial')
plt.plot(x, xs[num_steps], label='Final')

plt.xlabel('x')
plt.ylabel('y')

plt.legend()
plt.show()
