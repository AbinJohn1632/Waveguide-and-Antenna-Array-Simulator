#####TM#############

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# Parameters
a = 5  # Waveguide width
b = 2  # Waveguide height
c= 10 #depth
m, n = 4, 3  # TM_mn mode (change for other modes)
f = 10e9  # Frequency (Hz)
epsilon = 8.854e-12  # Permittivity
mu = 4e-7 * np.pi  # Permeability

lambda_ = (3 * (10**8)) / f
eta = np.sqrt(mu / epsilon)
omega = 2 * np.pi * f
k = omega * np.sqrt(mu * epsilon)

H=(np.pi * m / a)**2 +(np.pi * n / b)**2

# Characteristic equation and solve for beta
def char_eqn(beta):
    return k**2 -(np.pi * m / a)**2 - (np.pi * n / b)**2 -beta**2

beta = fsolve(char_eqn, 0)[0]

# Define field expressions for Ez
def Ez(x, y, z,t):
    return np.sin(m * np.pi * x / a) * np.sin(n * np.pi * y / b) * np.exp(1j*omega*t-1* beta * z)
def Ex(x, y, z, t):
    return (-beta/H)*np.cos(m * np.pi * x / a) * np.sin(n * np.pi * y / b) * np.exp(1j*omega*t-1* beta * z)
def Ey(x, y, z,t):
    return (-beta/H)*np.sin(m * np.pi * x / a) * np.cos(n * np.pi * y / b) * np.exp(1j*omega*t-1* beta * z)
def Hx(x, y, z,t):
    return (1j*omega*epsilon/H)*np.sin(m * np.pi * x / a) * np.cos(n * np.pi * y / b) * np.exp(1j*omega*t-1* beta * z)
def Hy(x, y, z,t):
    return (1j*omega*epsilon/H)*np.cos(m * np.pi * x / a) * np.sin(n * np.pi * y / b) * np.exp(1j*omega*t-1* beta * z)

# Create a meshgrid for the waveguide
x = np.linspace(0, a, 100)
y = np.linspace(0, b, 100)
z = np.linspace(0, c, 50)  # Depth from 0 to 1, you can adjust this based on your requirements
t = np.linspace(0, 100, 50)
X, Y, Z, T = np.meshgrid(x, y, z, t)

# Evaluate Hz at each point in the meshgrid
ez_values = Ez(X, Y, Z, T)
hx_values = Hx(X, Y, Z, T)
hy_values = Hy(X, Y, Z, T)
ex_values = Ex(X, Y, Z, T)
ey_values = Ey(X, Y, Z, T)

# Create a 3D plot
fig = plt.figure(figsize=(25,25))
ax = fig.add_subplot(241, projection='3d')

# Add quiver plot subplot
ax_quiver = fig.add_subplot(221)
ax_quiver.quiver(X[:,:,0,0], Y[:,:,0,0], hx_values[:,:,0,0], hy_values[:,:,0,0],  color='black')
ax_quiver.set_xlabel('X (m)')
ax_quiver.set_ylabel('Y (m)')
ax_quiver.set_xlim([0, a])
ax_quiver.set_ylim([0, b])
ax_quiver.set_title('Magnetic Field Vector(end view)')

# Add subplots for E feild
ax = fig.add_subplot(222)
ax.quiver(X[:,:,0,0], Y[:,:,0,0], ex_values[:,:,0,0], ey_values[:,:,0,0],color='black')  
ax.set_xlabel('Hx')
ax.set_ylabel('Hy')
ax.set_xlim([0, a])
ax.set_ylim([0, b])
ax.set_title('3D Quiver Plot ')

# 3D surface plot for Ex
ax_surf_ex = fig.add_subplot(223, projection='3d')
surf_ex = ax_surf_ex.plot_surface(X[:, :, 0, 0], Y[:, :, 0, 0], ez_values[:, :, 0, 0], cmap=cm.viridis, rstride=1, cstride=1, linewidth=0)
ax_surf_ex.set_xlabel('X (m)')
ax_surf_ex.set_ylabel('Y (m)')
ax_surf_ex.set_zlabel('Ex (V/m)')
ax_surf_ex.set_title('TM_{}_{} mode - Ex'.format(m, n))

ax_end = fig.add_subplot(224)
ax_end.contourf(X[:,:,0,0], Y[:,:,0,0], ez_values[:,:,0,0].real, cmap=cm.cividis)
ax_end.set_xlabel('X (m)')
ax_end.set_ylabel('Y (m)')
ax_end.set_title('End View of Magnetic Field (XY-plane)')

plt.tight_layout()
plt.show()