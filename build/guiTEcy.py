from scipy.optimize import fsolve,root
from matplotlib import cm
from scipy.special import jn, jn_zeros
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, Label,Scrollbar
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage, filedialog,ttk, Label, Scrollbar
import threading

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"E:\python\machine_learning\college_project\Waveguide_Simulator_sem6\build\assets\frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

window = Tk()
window.title("TE_circular")
window.geometry("1305x770")
window.configure(bg="#F1F3F4")

canvas = Canvas(
    window,
    bg="#F1F3F4",
    height=770,
    width=1305,
    bd=0,
    highlightthickness=0,
    relief="ridge"
)

canvas.place(x=0, y=0)
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

label_Head = tk.Label(window, text="Circular Waveguide(TE)", bg="#D9D9D9",font=(10), relief="sunken")
label_Head.place(x=0, y=20,width=370)

# Parameters
default_R = 0.01
default_f = 8.79e9
default_T = 1
default_Z= 1

c = 100  # Depth
m = 1 # TM_mn mode (m is an integer for cylindrical symmetry)
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

def char_eqn_scaled(beta):
    # Scale the equation if necessary
    return (beta**2 - k**2) * epsilon - k_c**2

result = root(char_eqn_scaled, x0=0.1 ,method='hybr')
beta = result.x[0]

A=1
k_g=beta

# Dropdown menu for selecting m and n values
m_values = list(range(6))

selected_m = tk.StringVar()
m_dropdown = ttk.Combobox(window, textvariable=selected_m, values=m_values)
m_dropdown.place(x=255, y=220, width=30)
label_m = tk.Label(window, text="m=", bg="#D9D9D9")
label_m.place(x=229, y=220)

# Define a function to update the equations when m and n change
def update_equations():
    global hr,h0,hz,er,e0,ez
    hr = Hr(R, Theta, Z,T)
    # hr=hr.flatten()
    h0 = Htheta(R, Theta, Z,T)
    # h0=h0.flatten()
    hz = Hz(R, Theta, Z,T)
    # hz=hz.flatten()
    er = Er(R, Theta, Z,T)
    # er=er.flatten()
    e0 = Etheta(R, Theta, Z,T)
    # e0=e0.flatten()
    ez = Ez(R, Theta, Z,T)
    # ez=ez.flatten()
    update_plot(plot_dropdown.get())

# Function to update m and n values
def update_m_n_values():
    global m,k_c
    m = int(selected_m.get())
    k_c = jn_zeros(m, 1)[0] / default_R

# Set default values for m
m_dropdown.set(str(m))

# Define field expressions for Ez
# Magnetic field components
def Hr(r, theta, z,t):
    return -(k_g * k_c /(omega*mu)) * A * jn(m, k_c * r) * np.cos(m * theta)* np.exp(1j*(-omega*t+ k_g * z))

def Htheta(r, theta, z,t):
    return (k_g / (omega *mu))* (m * A /r) * jn(m, k_c * r) * np.sin(m * theta)* np.exp(1j*(-omega*t+ k_g * z))

def Hz(r, theta, z,t):
    return 1j*(k_c**2 / (omega *mu))* A * jn(m, k_c * r) * np.cos(m * theta)* np.exp(1j*(-omega*t+k_g * z))  

# Electric field components
def Er(r, theta, z,t):
    return (m * A * jn(m, k_c * r) * np.sin(m * theta)/r)* np.exp(1j*(-omega*t+k_g * z))

def Etheta(r, theta, z,t):
    return k_c * A * jn(m, k_c * r) * np.cos(m * theta)* np.exp(1j*(-omega*t+k_g * z))

def Ez(r, theta, z,t):
    return np.zeros_like(r) 

window.update_equations = update_equations

label_T = tk.Label(window, text="Time Axis(t)", bg="#D9D9D9")
label_T.place(x=19, y=110)

entry_T = ttk.Scale(window, from_=1, to=19, orient="horizontal", length=320)
entry_T.place(x=19, y=130)
entry_T.set(default_T)

label_Z = tk.Label(window, text="Space Axis(z)", bg="#D9D9D9")
label_Z.place(x=19, y=160)

entry_Z = ttk.Scale(window, from_=1, to=19, orient="horizontal", length=320)
entry_Z.place(x=19, y=180)
entry_Z.set(default_Z)

# Dropdown menu for selecting plots
plot_options = ['Radial', 'Angular', 'End_View(z)']
selected_plot = tk.StringVar()
selected_plot.set('End_View(z)')
plot_dropdown = ttk.Combobox(window, textvariable=selected_plot, values=plot_options)
label_drop = tk.Label(window, text="Mode=", bg="#D9D9D9")
label_drop.place(x=19, y=220)
plot_dropdown.place(x=79, y=220)

# Generate polar grid points
Nr, Ntheta, Nz,Nt = 60,60,20,20
r = np.linspace(0.001, default_R,Nr)
theta = np.linspace(0, 2 * np.pi, Ntheta)
z = np.linspace(0, 100, Nz)
t = np.linspace(0, 100, Nt)
R, Theta, Z,T = np.meshgrid(r, theta, z,t)

def update_plot(selected_plot):
    plt.close('all')  # Close all previously opened figures

    T_new = int(entry_T.get())
    Z_new = int(entry_Z.get())

    fig = plt.figure(figsize=(9.1175, 5.3708333333))
    ax0=fig.add_subplot(121, polar=True)
    ax1=fig.add_subplot(122, polar=True)
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas_widget = canvas.get_tk_widget().place(x=379, y=14)

    if selected_plot == 'Radial':
        ax0.contourf(Theta[:, :, Z_new, T_new], R[:, :, Z_new, T_new], er[:, :, Z_new, T_new].real, cmap=cm.coolwarm)
        ax0.set_title('Radial Electric Field Component')
        ax0.grid(0)
        ax0.set_rticks([default_R]) 
        ax0.set_rmax(default_R)
        ax1.contourf(Theta[:, :, Z_new, T_new], R[:, :, Z_new, T_new], hr[:, :, Z_new, T_new].real, cmap=cm.coolwarm)
        ax1.set_title('Radial Magnetic Field Component')
        ax1.grid(0)
        ax1.set_rticks([default_R]) 
        ax1.set_rmax(default_R)

    elif selected_plot == 'Angular':
        ax0.contourf(Theta[:, :, Z_new, T_new], R[:, :, Z_new, T_new], e0[:, :, Z_new, T_new].real, cmap=cm.coolwarm)
        ax0.set_title('Angular Electric Field Component')
        ax0.grid(0)
        ax0.set_rticks([default_R]) 
        ax0.set_rmax(default_R)
        ax1.contourf(Theta[:, :, Z_new, T_new], R[:, :, Z_new, T_new], h0[:, :, Z_new, T_new].real, cmap=cm.coolwarm)
        ax1.set_title('Angular Electric Field Component')
        ax1.grid(0)
        ax1.set_rticks([default_R]) 
        ax1.set_rmax(default_R)

    elif selected_plot == 'End_View(z)':
        ax0.contourf(Theta[:, :, Z_new, T_new], R[:, :, Z_new, T_new], ez[:, :, Z_new, T_new].real, cmap=cm.coolwarm)
        ax0.set_title('Electric Field Component along Z axis')
        ax0.grid(0)
        ax0.set_rticks([default_R]) 
        ax0.set_rmax(default_R)
        ax1.contourf(Theta[:, :, Z_new, T_new], R[:, :, Z_new, T_new], hz[:, :, Z_new, T_new].real, cmap=cm.coolwarm)
        ax1.set_title('Electric Field Component along Z axis')
        ax1.grid(0)
        ax1.set_rticks([default_R]) 
        ax1.set_rmax(default_R)

    fig.canvas.draw()

def plot_with_loading():

    loading_label = tk.Label(window, text="Loading...", bg="#F1F3F4",font=(14))
    loading_label.place(x=380, y=540)
    window.update()
    threading.Thread(target=update_equations()).start()
    loading_label.destroy()

# Bind the update function to dropdown selection
m_dropdown.bind("<<ComboboxSelected>>", lambda event=None: update_m_n_values())
button_update = tk.Button(window, text="Plot", command=lambda value=None:plot_with_loading())
button_update.place(x=125, y=260, width=100,height=50)
plot_dropdown.bind("<<ComboboxSelected>>", lambda event=None: update_plot(plot_dropdown.get()))
entry_T['command'] = lambda value=None: update_plot(plot_dropdown.get())
entry_Z['command'] = lambda value=None: update_plot(plot_dropdown.get())

import subprocess

def run_program_waveguide_TE():
    window.destroy()
    subprocess.call(["python", "build\guiTErec.py"])

# Button to switch to window 2
button_window2 = tk.Button(window, text="TE Rectangular", command=run_program_waveguide_TE)
button_window2.place(x=19, y=58, width=100,height=50)


def run_program_waveguide_TM():
    window.destroy()
    subprocess.call(["python", "build\guiTMrec.py"])

# Button to switch to window 2
button_to_TM = tk.Button(window, text="TM Rectangular", command=run_program_waveguide_TM)
button_to_TM.place(x=122, y=58, width=100,height=50)

def run_program_waveguide_TMcy():
    window.destroy()
    subprocess.call(["python", "build\guiTMcy.py"])

button_to_TMcy = tk.Button(window, text="TM Circular", command=run_program_waveguide_TMcy)
button_to_TMcy.place(x=225, y=58, width=100,height=50)

def to_home():
    window.destroy()
    subprocess.call(["python", "build\HOME.py"])

# Button to switch to window 2
button_To_home = tk.Button(window, text="HOME", command=to_home)
button_To_home.place(x=19, y=700, width=100,height=50)

def on_closing():
    window.destroy()  # Destroy the Tkinter window when closed
    exit()  # Exit the Python script, closing the terminal

window.protocol("WM_DELETE_WINDOW", on_closing)

window.resizable(0, 0)
window.mainloop()