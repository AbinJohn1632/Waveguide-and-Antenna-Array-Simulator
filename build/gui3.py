from scipy.optimize import fsolve,root
from matplotlib import cm
from scipy.special import jn, jn_zeros
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, Label,Scrollbar
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage,filedialog
from PIL import Image, ImageTk 

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"E:\python\machine_learning\college_project\Waveguide_Simulator_sem6\build\assets\frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

window = Tk()
window.title("TM_cylindrical")
window.geometry("1305x850")
window.configure(bg = "#F1F3F4")

canvas = Canvas(
    window,
    bg = "#F1F3F4",
    height = 850,
    width = 1305,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge"
)

canvas.place(x = 0, y = 0)
canvas.create_rectangle(
    0.0,
    17.0,
    370.0,
    831.0,
    fill="#D9D9D9",
    outline="")

canvas.create_rectangle(
    379.402015341185,
    568.0,
    1286.4012155902456,
    835.9999857976493,
    fill="#D9D9D9",
    outline="")

# Parameters
default_R = 0.01
default_f = 8.79e9
default_T = 0
default_Z= 10

c = 100  # Depth
m = 1 # TE_mn mode (m is an integer for cylindrical symmetry)
epsilon = 1  # Relative permittivity
mu = 4e-7 * np.pi  # Permeability

# Calculate angular frequency and wave number
omega = 2 * np.pi * default_f
k = omega * np.sqrt(mu*epsilon)

# Calculate the waveguide mode parameter kc
k_c = jn_zeros(m, 1)[0] / default_R

# Calculate angular frequency and wave number for default parameters
omega_default = 2 * np.pi * default_f
k_default = omega_default * np.sqrt(mu * epsilon)

# Calculate the waveguide mode parameter kc for default parameters
k_c_default = jn_zeros(m, 1)[0] / default_R

# def char_eqn(beta):
#     return (beta**2 - k_default**2) * epsilon_r - k_c_default**2
# beta = fsolve(char_eqn, 0)[0]

def char_eqn_scaled(beta):
    # Scale the equation if necessary
    return (beta**2 - k**2) * epsilon - k_c**2

result = root(char_eqn_scaled, x0=0.1 ,method='hybr')
beta = result.x[0]

A=1
k_g=beta

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

fig = plt.figure(figsize=(9.1175, 5.3708333333))

ax0 = fig.add_subplot(231, polar=True)
ax1 = fig.add_subplot(232, polar=True)
ax2 = fig.add_subplot(233, polar=True)
ax3 = fig.add_subplot(234, polar=True)
ax4 = fig.add_subplot(235, polar=True)
ax5 = fig.add_subplot(236, polar=True)

canvas = FigureCanvasTkAgg(fig, master=window)
canvas_widget = canvas.get_tk_widget().place(x=379,y=14)

label_T = tk.Label(window, text="T", bg="#D9D9D9")
label_T.place(x=19, y=160)

entry_T = ttk.Scale(window, from_=1, to=50, orient="horizontal", length=320)
entry_T.place(x=19, y=180)
entry_T.set(default_T)

label_Z = tk.Label(window, text="Z", bg="#D9D9D9")
label_Z.place(x=19, y=110)

entry_Z = ttk.Scale(window, from_=1, to=50, orient="horizontal", length=320)
entry_Z.place(x=19, y=130)
entry_Z.set(default_Z)

# Generate polar grid points
Nr, Ntheta, Nz,Nt = 50, 50, 50,50
r = np.linspace(0.0001, 0.01, Nr)
theta = np.linspace(0, 2 * np.pi, Ntheta)
z = np.linspace(0, 100, Nz)
t = np.linspace(0, 100, Nt)
R, Theta, Z,T = np.meshgrid(r, theta, z,t, indexing='ij')

# Calculate field magnitudes at grid points
er = Hr(R, Theta, Z,T)
etheta = Htheta(R, Theta, Z,T)
ez = Hz(R, Theta, Z,T)
hr = Er(R, Theta, Z,T)
htheta = Etheta(R, Theta, Z,T)
hz = Ez(R, Theta, Z,T)

def update_plot(event=None):

    T_new = int(entry_T.get())
    Z_new = int(entry_Z.get())

    ax0.clear()
    ax0.contourf(Theta[:, :, Z_new,T_new], R[:, :, Z_new,T_new], er[:, :, Z_new,T_new].real, cmap=cm.coolwarm)
    ax0.set_title('Hr')

    ax1.clear()
    ax1.contourf(Theta[:, :, Z_new,T_new], R[:, :, Z_new,T_new], etheta[:, :, Z_new,T_new].real, cmap=cm.coolwarm)
    ax1.set_title('Htheta')

    ax2.clear()
    ax2.contourf(Theta[:, :, Z_new,T_new], R[:, :, Z_new,T_new], ez[:, :, Z_new,T_new].real, cmap=cm.coolwarm)
    ax2.set_title('Hz')

    ax3.clear()
    ax3.contourf(Theta[:, :, Z_new,T_new], R[:, :, Z_new,T_new], hr[:, :, Z_new,T_new].real, cmap=cm.coolwarm)
    ax3.set_title('Er')

    ax4.clear()
    ax4.contourf(Theta[:, :, Z_new,T_new], R[:, :, Z_new,T_new], htheta[:, :, Z_new,T_new].real, cmap=cm.coolwarm)
    ax4.set_title('Etheta')

    ax5.clear()
    ax5.contourf(Theta[:, :, Z_new,T_new], R[:, :, Z_new,T_new], hz[:, :, Z_new,T_new].real, cmap=cm.coolwarm)
    ax5.set_title('Ez')

    canvas.draw()


entry_T['command']=update_plot
entry_Z['command']=update_plot

window.resizable(0,0)
window.mainloop()



