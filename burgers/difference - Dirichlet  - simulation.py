#---------------------------------------------------------
# This program uses finitie difference method to simulate
# Burgers equation.
#---------------------------------------------------------
# By Bowen Li
# October 7, 2024
#---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
import time


# The iteration map x_{i} = f(x_{i-1})
def f(x, dt):
   """brugers model"""
   # Setting up vector
   d = np.zeros(N)
   # Loops over indices (with operations and Python underflow indexing handling edge cases)
   d[0] = 0
   d[N - 1] = 0
   for i in range(1, N-1):
      d[i] = 0.5*x[i+1] + 0.5*x[i-1] - 0.25*(dt/dx)*(x[i+1]**2 - x[i-1]**2) + nu*(dt/(dx**2))*(x[i+1] - 2*x[i] + x[i-1])
   return d


time_start = time.time()

# These are our constants
N = 101  # Number of variables

nu = 0.005

dt = 0.02
num_steps = 100
dx = np.pi/(N-1)

x = np.linspace(0, np.pi, N)
x0 = np.sin(x)  # Initial state (equilibrium)
#x0 += np.random.normal(0, 1, N)  # Add small perturbation to the first variable

'''t_span = [0.0, num_steps * dt]
#t = np.arange(0.0, (num_steps + 1)*dt, dt)
t = np.linspace(0.0, num_steps * dt, num_steps + 1)
x = np.linspace(-np.pi, np.pi, num_steps + 1)

xs = solve_ivp(g, t_span, x0, t_eval=t).y
xs = xs.T'''

# Runge-Kutta method
xs = np.zeros([num_steps+1, N])
xs[0] = x0
print(xs[0])
for i in range(num_steps):
   xs[i+1] = f(xs[i], dt)


plt.figure(1)
plt.plot(x, xs[0], label='Initial')
plt.plot(x, xs[num_steps], label='Final')

plt.xlabel('x')
plt.ylabel('y')

plt.legend()
plt.show()
