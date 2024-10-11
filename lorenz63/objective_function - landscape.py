#-----------------------------------------------------------------
# This program draws the figure of the objective function at z=20.
#-----------------------------------------------------------------
# By Bowen Li
# October 7, 2024
#-----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import time


def g(x, t):

   """lorenz 63 model"""
   # Setting up vector
   d = np.zeros(3)
   # Loops over indices (with operations and Python underflow indexing handling edge cases)
   d[0] = ss*(x[1] - x[0])
   d[1] = r*x[0] - x[1] - x[0]*x[2]
   d[2] = x[0]*x[1] - b*x[2]
   return d

# The iteration map x_{i} = f(x_{i-1})
def f(x, dt):
   t = 0
   k1 = g(x, t)
   k2 = g(x + 0.5*dt*k1, t)
   k3 = g(x + 0.5*dt*k2, t)
   k4 = g(x + dt*k3, t)
   return x + dt*(k1 + 2*k2 + 2*k3 + k4)/6


def obj(x0):
   fv = 0.5*alpha*np.linalg.norm(x0 - xs[0])**2
   xt = odeint(g, x0, t)
   for k in range(num_steps//M + 1):
      fv += (0.5/(num_steps//M)) * np.linalg.norm(xt[k*M] - xs[k*M])**2
   return fv


# These are our constants
ss=10
r=28
b=8/3


x0 = [-0.5, 0.5, 20.5]  # Initial state (equilibrium)
#x0[0] += 0.01 #np.random.normal(0, 1, N)  # Add small perturbation to the first variable

dt = 0.01
num_steps = 300
alpha = 0.1

t = np.arange(0.0, (num_steps + 1)*dt, dt)

xs = odeint(g, x0, t)
print(xs[0])

M = num_steps//10


time_start = time.time()

H = 100

xx = np.linspace(-5,5,H)
yy = np.linspace(-5,5,H)
X, Y = np.meshgrid(xx, yy)
print(X.shape)

obj_f = np.zeros([H, H])
for i in range(H):
   for j in range(H):
      obj_f[i, j] = obj([X[i,j], Y[i,j], 20])

fig = plt.figure(dpi=60, figsize=(15, 15))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, obj_f, cmap='viridis', shade=True)

# Adjust the perspective
ax.view_init(elev=30, azim=320)  # elev is the elevation angle, and azim is the azimuth angle

# Add a color bar
mappable = plt.cm.ScalarMappable(cmap='viridis')
mappable.set_array([0,400])
cbar = fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=10, pad=0.1)

#ax.plot(x0[0], x0[1], 0, c='r', marker='x')
#ax.plot(0, 0, obj(np.zeros(3)), c='g', marker='x')

time_end = time.time()
print("Time: "+str(time_end - time_start)+" Seconds")

# Adjust the distance of the X, Y, and Z axis tick labels
ax.xaxis.set_tick_params(pad=6)  # Adjust the distance between the X-axis tick labels and the axis
ax.yaxis.set_tick_params(pad=6)  # Adjust the distance between the Y-axis tick labels and the axis
ax.zaxis.set_tick_params(pad=10)  # Adjust the distance between the Z-axis tick labels and the axis

ax.set_zticks(np.linspace(0, 400, 5))

# Manually set the ticks for the color bar
cbar.set_ticks(np.linspace(0, 400, 5))
cbar.update_ticks()

# Set the font size of the color bar ticks
cbar.ax.tick_params(labelsize=21)

#ax.set_title('3D Surface plot', fontsize=18)
ax.set_xlabel('X', fontsize=25, labelpad=16)
ax.set_ylabel('Y', fontsize=25, labelpad=16)
ax.set_zlabel('Z', fontsize=25, labelpad=16)

ax.tick_params(axis='both', which='major', labelsize=21)
ax.tick_params(axis='z', labelsize=21)

# Set the background to white
ax.set_facecolor('white')

# Remove the grid lines
ax.grid(False)

plt.savefig('landscape.png', dpi=600, format='png', bbox_inches='tight')  # Use high DPI when saving, then use an image tool to convert the PNG to PDF to reduce the size of figure

#plt.show()
