#---------------------------------------------------------
# This program verifies the tangent linear model.
#---------------------------------------------------------
# By Bowen Li
# October 7, 2024
#---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import odeint

'''
# Spectral
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

# The Jacobi of g(x)
def dg(x):
   A = np.zeros((N, N))
   for i in range(N):
      for n in range(i):
         A[i, n] -= 0.5*(i - n)*x[i - n - 1]
         A[i, i - n - 1] -= 0.5*(i - n)*x[n]
      for n in range(i+1, N):
         A[i, n] -= 0.5*(n - i)*x[n - i - 1]
         A[i, n - i - 1] -= 0.5*(n - i)*x[n]
      for n in range(N-i-1):
         A[i, n] += 0.5*(n + i + 2)*x[n + i + 1]
         A[i, n + i + 1] += 0.5*(n + i + 2)*x[n]
      A[i, i] -= nu * ((i+1)**2)
   return A

'''
# FEM
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

# The Jacobi of g(x)
def dg(x):
   A = np.zeros((N, N))
   #x[0] = 0
   #x[N - 1] = 0
   for i in range(1, N-1):
      A[i, i-1] = (1/6) * (2*x[i-1] + x[i]) + (nu/(dx))
      A[i, i] = (1/6) * (x[i-1] - x[i+1]) - (2*nu/(dx))
      A[i, i+1] = (1/6) * (-x[i] - 2*x[i+1]) + (nu/(dx))
   return np.matmul(G, A)

'''
# FD
def g(x, t):

   """brugers model"""
   # Setting up vector
   d = np.zeros(N)
   # Loops over indices (with operations and Python underflow indexing handling edge cases)
   d[0] = 0
   d[N - 1] = 0
   for i in range(1, N-1):
      d[i] = 0.5*x[i+1] + 0.5*x[i-1] - 0.25*(dt/dx)*(x[i+1]**2 - x[i-1]**2) + nu*(dt/(dx**2))*(x[i+1] - 2*x[i] + x[i-1])
   return d

# Tangent linear model
def dg(x):
   A = np.zeros([N, N])
   for i in range(1, N-1):
      A[i, i-1] = 0.5 + 0.5*(dt/dx)*x[i-1] + nu*(dt/(dx**2))
      A[i, i] = -2*nu*(dt/(dx**2))
      A[i, i+1] = 0.5 - 0.5*(dt/dx)*x[i+1] + nu*(dt/(dx**2))
   return A
'''

N = 100  # Number of variables


nu = 0.005
x0 = np.zeros(N)  # Initial state (equilibrium)
#x0 += np.random.normal(0, 1, N)  # Add small perturbation to the first variable
x0[0] = 1

dt = 0.01
num_steps = 200
alpha = 0.1

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

t = np.arange(0.0, (num_steps + 1)*dt, dt)

x = np.random.randn(N)

dd = np.random.randn(N)
dd = dd/np.linalg.norm(dd)

# Calculate the relative err
eps = 1e-7
approx = (g(x + eps*dd, t) - g(x - eps*dd, t))/(2*eps)
real = np.matmul(dg(x), dd)

err = np.linalg.norm(approx - real)/(np.linalg.norm(approx) + np.linalg.norm(real))

print(err)