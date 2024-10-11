#-----------------------------------------------
# This program uses spectral method to simulate
# Burgers equation.
#-----------------------------------------------
# By Bowen Li
# October 7, 2024
#-----------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
import time


#def g(x, t):
#
#   """burgers model"""
#   # Setting up vector
#   d = np.zeros(N)
#   # Loops over indices (with operations and Python underflow #indexing handling edge cases)
#   for i in range(N):
#      for n in range(N):
#         m = i - n
#         if m >= 1 and m <= N:
#            d[i] -= m * x[n] * x[m-1] / 2
#         m = i + n + 2 
#         if m >=1 and m <= N:
#            d[i] += m * x[n] * x[m-1] / 2
#         m = n - i
#         if m >=1 and m <= N:
#            d[i] -= m * x[n] * x[m-1] / 2
#      d[i] -= nu * ((i+1)**2) * x[i]
#   return d

def g(x, t):

   """brugers model"""
   # Setting up vector
   d = np.zeros(N)
   # Loops over indices (with operations and Python underflow indexing handling edge cases)
   for i in range(N):
      for n in range(i):
         d[i] -= 0.5*(i - n)*x[n]*x[i - n - 1]
      for n in range(i+1, N):
         d[i] -= 0.5*(n - i)*x[n]*x[n - i - 1]
      for n in range(N-i-1):
         d[i] += 0.5*(n + i + 2)*x[n]*x[n + i + 1]
      d[i] -= nu * ((i+1)**2) * x[i]
   return d

# The iteration map x_{i} = f(x_{i-1})
def f(x, dt):
   t = 0
   return x + dt*g(x, t)


time_start = time.time()

# These are our constants
N = 100  # Number of variables

nu = 0.005
x0 = np.zeros(N)  # Initial state (equilibrium)
#x0 += np.random.normal(0, 1, N)  # Add small perturbation to the first variable
x0[0] = 1

dt = 0.01
num_steps = 200

'''t_span = [0.0, num_steps * dt]
#t = np.arange(0.0, (num_steps + 1)*dt, dt)
t = np.linspace(0.0, num_steps * dt, num_steps + 1)
x = np.linspace(-np.pi, np.pi, num_steps + 1)

xs = solve_ivp(g, t_span, x0, t_eval=t).y
xs = xs.T'''

# Time integration
x = np.linspace(0, np.pi, 100)
xs = np.zeros([num_steps+1, N])
xs[0] = x0
for i in range(num_steps):
   xs[i+1] = f(xs[i], dt)

print(xs[0])

FM = np.zeros([N, 100])
for i in range(N):
   FM[i] = np.sin((i+1)*x)
FM = FM.T


plt.figure(1)
plt.plot(x, np.matmul(FM, xs[0]), label='Initial')
plt.plot(x, np.matmul(FM, xs[num_steps]), label='Final')

plt.xlabel('x')
plt.ylabel('y')

plt.legend()
plt.show()
