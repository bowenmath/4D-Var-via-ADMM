#---------------------------------------------------------
# This program uses linearized multi-block ADMM to recover 
# Burgers equation with finitie element discretization.
#---------------------------------------------------------
# By Bowen Li
# October 7, 2024
#---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import odeint
import time


def g(x, t):

   """brugers model with finite element method"""
   # Setting up vector
   d = np.zeros(N)
   # Loops over indices (with operations and Python underflow indexing handling edge cases)
   d[0] = 0
   d[N - 1] = 0
   for i in range(1, N-1):
      d[i] = (1/6) * (x[i-1]**2 + x[i]*x[i-1] - x[i+1]*x[i] - x[i+1]**2) + (nu/(dx)) * (x[i-1] - 2*x[i] + x[i+1])
   d = np.matmul(G, d)
   d[0] = 0
   d[N - 1] = 0
   return d

# The Jacobi of g(x)
def dg(x):
   A = np.zeros((N, N))
   for i in range(1, N-1):
      A[i, i-1] = (1/6) * (2*x[i-1] + x[i]) + (nu/(dx))
      A[i, i] = (1/6) * (x[i-1] - x[i+1]) - (2*nu/(dx))
      A[i, i+1] = (1/6) * (-x[i] - 2*x[i+1]) + (nu/(dx))
   return np.matmul(G, A)

# The iteration map x_{i} = f(x_{i-1})
def f(x, dt):
   t = 0
   return x + dt*g(x, t)

# Tangent linear model
def df(x, dt):
   return np.eye(N) + dt*dg(x)

# Solve the subproblems of ADMM
def subproblem(x0, x1, x2, lmda1, lmda2, xs, dt, rho, *, l=1, m=1, n=1):
   # The linearized subproblem
   # Use l,m,n to control which term appears in the subproblem
   # x1 is the variable to update

   eta = 0.1  # The parameter \eta for linearization in ADMM
   para = 1 + eta * (2*l + rho * m)
   b = 2*l*xs + m * (rho * f(x0, dt) + lmda1) + n * np.matmul(df(x1, dt).T, rho * x2 - rho*f(x1, dt) - lmda2)
   x1 = (x1 + eta * b)/para
   return x1

# Evaluate the objective function given the initial value
def obj(x0):
   fv = 0.5*alpha*np.linalg.norm(x0 - xs[0])**2
   xt = odeint(g, x0, t)
   for k in range(num_steps//M + 1):
      fv += (0.5/(num_steps//M)) * np.linalg.norm(xt[k*M] - xs[k*M])**2
   return fv


time_start = time.time()  # Record time

# These are our constants
N = 101  # Number of variables


nu = 0.005  # The viscosity coefficient

x = np.linspace(0, np.pi, N)
x0 = np.sin(x)  # Initial state (equilibrium)

dt = 0.01  # The size of time steps
num_steps = 200  # The number of time steps
T = dt*num_steps  # Total time length
alpha = 0.1  # The parameter for background error (Here we set the background information equal to the inital observation for simplicity)

dx = np.pi/(N-1)  # The size of space discretization


# The simulation of burgers model

# Generate the stiffness matrix
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


# Generate the matrix for initial value
B = np.zeros([N, N])
for i in range(1, N-1):
   B[i,i-1] = -1/dx
   B[i, i] = 2/dx
   B[i, i+1] = -1/dx
x0 = np.matmul(G, np.matmul(B, x0))
x0[0] = 0
x0[N-1] = 0


t = np.arange(0.0, (num_steps + 1)*dt, dt)

#xst = odeint(g, x0, t)
xst = np.zeros([num_steps + 1, N])
xst[0] = x0
for i in range(num_steps):
   xst[i + 1] = f(xst[i], dt)

print(xst[0])


iter_steps = 1000  # Number of iterations for ADMM
M = num_steps//10  # Take observations ever M steps (Totally 10 observations)
rho = 1.5  # The parameter \rho in ADMM
mu = 20  # The scaling parameter \mu in ADMM


x = np.zeros((num_steps + 1, N))  # The main variable, corresponding to the variable u in the paper
xt = np.zeros((num_steps + 1, N))  # Intermediate variable for x
lmda = np.zeros((num_steps + 1, N))  # The array for \lambda (The first element is not used)
lmdat = np.zeros((num_steps + 1, N))  # Intermediate variable for \lambda
fval = np.zeros(iter_steps + 1)  # Store the value of objective function for each iteration
cons_err = np.zeros(iter_steps + 1)  # Store the constraint error for each iteration


# Add noise to the observations
seed = 0
np.random.seed(seed)

xs = np.zeros([num_steps+1, N])
for i in range(num_steps+1):
   if i % M == 0:
      xs[i] = xst[i] + 0.1*np.random.normal(0, 1, N)
print(xs[0])


# Calculate the value of objective function for initial value
fval[0] = dx*T*0.5*alpha*np.linalg.norm(x[0] - xs[0])**2
for k in range(num_steps//M + 1):
   fval[0] += (dx*T*0.5/(num_steps//M)) * np.linalg.norm(x[k*M] - xs[k*M])**2

# Calculate the constraint error for initial value
for k in range(num_steps):
   cons_err[0] += np.linalg.norm(x[k+1] - f(x[k], dt))**2


# The multi-block ADMM with Jacobian decomposition
for i in range(iter_steps):
   xt[0] = subproblem(x[0], x[0], x[1], lmda[1], lmda[1], xs[0], dt, rho, l=(0.5/(num_steps//M))*mu + 0.5*alpha*mu, m=0)
   for j in range(1, num_steps):
      if j % M == 0:
         xt[j] = subproblem(x[j-1], x[j], x[j+1], lmda[j], lmda[j+1], xs[j], dt, rho, l=(0.5/(num_steps//M))*mu)
      else:
         xt[j] = subproblem(x[j-1], x[j], x[j+1], lmda[j], lmda[j+1], xs[j], dt, rho, l=0)
      # Update \lambda
      lmdat[j] = lmda[j] - rho*(xt[j] - f(xt[j - 1], dt))

   xt[num_steps] = subproblem(x[num_steps-1], x[num_steps], x[num_steps], lmda[num_steps], lmda[num_steps], xs[num_steps], dt, rho, l=(0.5/(num_steps//M))*mu, n=0)
   # Update \lambda
   lmdat[num_steps] = lmda[num_steps] - rho*(xt[num_steps] - f(xt[num_steps - 1], dt))

   # Update all the variables to implement Jacobi decomposed ADMM
   for j in range(num_steps + 1):
      x[j] = xt[j]
      lmda[j] = lmdat[j]

   # Store the value of objective function
   fval[i + 1] = dx*T*0.5*alpha*np.linalg.norm(x[0] - xs[0])**2
   for k in range(num_steps//M + 1):
      fval[i + 1] += (dx*T*0.5/(num_steps//M)) * np.linalg.norm(x[k*M] - xs[k*M])**2
   # Store the constraint error
   for k in range(num_steps):
      cons_err[i + 1] += np.linalg.norm(x[k+1] - f(x[k], dt))**2


# Display the results
print("\n")
print("ADMM: ")
print(x[0])
print("Error:", fval[iter_steps])

print("Constraint Error:", cons_err[iter_steps])

time_end = time.time()
print("Time: "+str(time_end - time_start))


# Save an array as a text file
#np.savetxt('FEM_fv.txt', fval)


# Convergence graph of function value

xx = np.linspace(0, np.pi, N)
#xxt = np.zeros([num_steps+1, N])
#xxt[0] = x[0]
#for i in range(num_steps):
#   xxt[i+1] = f(xxt[i], dt)


plt.figure(0, dpi=60, figsize=(15, 10))
plt.plot(range(iter_steps + 1), fval, label='ADMM', color='#ff7f0e', linewidth=2)

# Adjust the font size of the tick marks
plt.tick_params(axis='both', which='major', labelsize=21)

#plt.axvline(x=200, color='k', linestyle='--')

plt.yscale('log')

plt.xlabel('Iteration', fontsize=26, labelpad=15)
plt.ylabel(r"$F - F^\star$", fontsize=26, labelpad=15)

plt.legend(fontsize=26, loc='upper right', frameon=True, framealpha=0.8, edgecolor='black')

plt.savefig('convergence_FEM.pdf', dpi=60, format='pdf', bbox_inches='tight')





# Convergence graph of constraint error
plt.figure(1, dpi=60, figsize=(15, 10))
plt.plot(range(1, iter_steps + 1), cons_err[1:], label='ADMM', color='#ff7f0e', linewidth=2)

# Adjust the font size of the tick marks
plt.tick_params(axis='both', which='major', labelsize=21)

#plt.axvline(x=200, color='k', linestyle='--')

plt.yscale('log')

plt.xlabel('Iteration', fontsize=26, labelpad=15)
plt.ylabel(r"$\sum_{k=0}^{N-1}\ \left\| \mathbf{u}_{k+1}-H(\mathbf{u}_k) \right\|^2$", fontsize=26, labelpad=15)

plt.legend(fontsize=26, loc='upper right', frameon=True, framealpha=0.8, edgecolor='black')

plt.savefig('cons_err_FEM.pdf', dpi=60, format='pdf', bbox_inches='tight')






# Create a 1x5 grid of subplot
fig, axs = plt.subplots(1, 4, figsize=(16, 4))

# Add data to each subplot
for i in range(4):
    axs[i].scatter(xx, xs[i*int(0.6/dt)], label='FE-Obs', c='magenta', marker='o', s=12)
    axs[i].plot(xx, xst[i*int(0.6/dt)], label='FE-True', color='#1f77b4', linewidth=3)
    axs[i].plot(xx, x[i*int(0.6/dt)], label='ADMM', linestyle='-.', color='#ff7f0e', linewidth=2)

    axs[i].set_title(f't = {round(i*0.6, 6)}', fontsize=14)
    axs[i].set_xlabel(f'x', fontsize=12)  # Set the x-axis label
    axs[i].set_ylabel(r"$\mathbf{u}$", fontsize=12)  # Set the y-axis label

    # Get the current subplot's handle and label
    handles, labels = axs[i].get_legend_handles_labels()

    # Manually adjust the order of handles and labels (e.g., if you need to change the order, you can reorder handles and labels)
    axs[i].legend([handles[0], handles[1], handles[2]], [labels[0], labels[1], labels[2]], fontsize=10, loc='upper left', frameon=True, framealpha=0.8, edgecolor='black')  # Show the legend

    # Set the y-axis range
    axs[i].set_ylim(-0.05, 1.55)

    axs[i].set_yticks(np.linspace(0.0, 1.5, 6))

    # Ensure each subplot has its own tick font size set individually
    axs[i].tick_params(axis='both', which='major', labelsize=10)  # Adjust the tick font size


plt.tight_layout()

# Save a high-quality image
plt.savefig('dynamic_FEM.pdf', dpi=60, format='pdf', bbox_inches='tight')

plt.show()