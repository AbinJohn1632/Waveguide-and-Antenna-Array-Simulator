import numpy as np
import matplotlib.pyplot as plt

def para(d, w,de):
    """Return the power as a function of azimuthal angle, phi."""
    phi = np.linspace(0, 2*np.pi, 1000)
    psi = (2*np.pi * d / lam) * np.cos(phi)+de
    return psi,phi

def normal_rad_pattern(psi,N):
    pattern=(np.sin(N*psi/2)/np.sin(psi/2))/N
    return pattern

# Wavelength, antenna spacing, feed coefficients.
de=0.1 #electric phase
lam = 1
d = lam*8
N = 2

# Calculate gain and directive gain; plot on a polar chart.
psi,phi = para(d, N,de)
pattern=normal_rad_pattern(psi,N)

fig = plt.figure(figsize=(10, 6))

ax=fig.add_subplot(121,polar=True)
ax.plot(phi, pattern, linewidth=1)
ax.set_title("Polar plot Uniform linear array for N={}".format(N))

ax2=fig.add_subplot(122)
ax2.plot(phi, pattern, linewidth=1)
ax2.set_xlabel("phi")
ax2.set_ylabel("pattern")
ax2.set_title("Uniform linear array for N={}".format(N))


plt.tight_layout()
plt.show()
