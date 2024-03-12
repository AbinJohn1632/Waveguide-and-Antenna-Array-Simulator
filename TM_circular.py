import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,root
from matplotlib import cm
from scipy.special import jn, jn_zeros

# Parameters
R = 0.01  # Inner radius of the waveguide (modify as needed)
c = 10  # Depth
m = 1 # TM_mn mode (m is an integer for cylindrical symmetry)
f = 8.79e9  # Cut-off frequency for TE11 mode (Hz) - modify as needed
epsilon = 1  # Relative permittivity
mu = 4e-7 * np.pi  # Permeability

# Calculate angular frequency and wave number
omega = 2 * np.pi * f
k = omega * np.sqrt(mu*epsilon)

# Calculate the waveguide mode parameter kc
k_c = jn_zeros(m, 1)[0] / R

def char_eqn_scaled(beta):
    # Scale the equation if necessary
    return (beta**2 - k**2) * epsilon - k_c**2

# beta = fsolve(char_eqn_scaled, 0.1)[0]
result = root(char_eqn_scaled, x0=0.1 ,method='hybr')
beta = result.x[0]

k_g = beta

A = 1

# Define field expressions for Az in cylindrical coordinates
def Az(r, theta, z,t):
    return A * jn(m, k_c * r) * np.cos(m * theta) * np.exp(1j *( beta * z-omega*t))

# Magnetic field components
def Hr(r, theta, z,t):
    return -m * A * (jn(m, k_c * r)/r) * np.sin(m * theta) * np.exp(1j *( beta * z-omega*t))

def Htheta(r, theta, z,t):
    return -k_c * A * jn(m, k_c * r) * np.cos(m * theta) * np.exp(1j *( beta * z-omega*t))

def Hz(r, theta, z,t):
    return np.zeros_like(r)  # Hz is zero in the provided equations

# Electric field components
def Er(r, theta, z,t):
    return -k_g  * (k_c/(epsilon * omega)) * A * jn(m, k_c * r) * np.cos(m * theta) * np.exp(1j *( beta * z-omega*t))

def Etheta(r, theta, z,t):
    return (k_g /(epsilon*omega*r)) * m * A * jn(m, k_c * r) * np.sin(m * theta) * np.exp(1j *( beta * z-omega*t))

def Ez(r, theta, z,t):
    return (-1j * A /(r**2 * epsilon * omega))* (k**2 * r**2 * jn(m, k_c * r).real + k_c * r * jn(m, k_c * r).imag - m**2 * jn(m, k_c * r).real) * np.cos(m * theta) * np.exp(1j *( beta * z-omega*t))

# Generate polar grid points
Nr, Ntheta, Nz,Nt = 50, 50, 50,50
r = np.linspace(0.00001, R, Nr)
theta = np.linspace(0, 2 * np.pi, Ntheta)
z = np.linspace(0, c, Nz)
t = np.linspace(0, 100, Nt)
R, Theta, Z,T = np.meshgrid(r, theta, z,t, indexing='ij')

# Calculate field magnitudes at grid points
er = Hr(R, Theta, Z,T)
etheta = Htheta(R, Theta, Z,T)
ez = Hz(R, Theta, Z,T)
hr = Er(R, Theta, Z,T)
htheta = Etheta(R, Theta, Z,T)
hz = Ez(R, Theta, Z,T)

fig = plt.figure(figsize=(9.1175, 5.3708333333))

ax0 = fig.add_subplot(231, polar=True)
ax1 = fig.add_subplot(232, polar=True)
ax2 = fig.add_subplot(233, polar=True)
ax3 = fig.add_subplot(234, polar=True)
ax4 = fig.add_subplot(235, polar=True)
ax5 = fig.add_subplot(236, polar=True)

# ax0.clear()
ax0.contourf(Theta[:, :, 0,33], R[:, :, 0,33], er[:, :, 0,33].real, cmap=cm.viridis)
ax0.set_title('Hr')

# ax1.clear()
ax1.contourf(Theta[:, :, 0,33], R[:, :, 0,33], etheta[:, :, 0,33].real, cmap=cm.viridis)
ax1.set_title('Htheta')

# ax2.clear()
ax2.contourf(Theta[:, :, 0,33], R[:, :, 0,33], ez[:, :, 0,33].real, cmap=cm.viridis)
ax2.set_title('Hz')

# ax3.clear()
ax3.contourf(Theta[:, :, 0,33], R[:, :, 0,33], hr[:, :, 0,33].real, cmap=cm.viridis)
ax3.set_title('Er')

# ax4.clear()
ax4.contourf(Theta[:, :, 0,33], R[:, :, 0,33], htheta[:, :, 0,33].real, cmap=cm.viridis)
ax4.set_title('Etheta')

# ax5.clear()
ax5.contourf(Theta[:, :, 0,33], R[:, :, 0,33], hz[:, :, 0,33].real, cmap=cm.viridis)
ax5.set_title('Ez')

plt.tight_layout()
plt.show()
