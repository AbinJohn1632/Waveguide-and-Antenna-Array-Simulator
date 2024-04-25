from scipy.optimize import root
from matplotlib import cm
from scipy.special import jn, jn_zeros
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage, filedialog,ttk, Label, Scrollbar
import threading

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"E:\python\machine_learning\college_project\Waveguide_Simulator_sem6\build\assets\frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

window = Tk()
window.title("TM_rectangular")
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

label_Head = tk.Label(window, text="Rectangular Waveguide(TM)", bg="#D9D9D9",font=(10), relief="sunken")
label_Head.place(x=0, y=20,width=370)

# Parameters
default_f = 8.79e9
default_T = 1
default_Z = 1

# Parameters
a = 5  # Waveguide width
b = 2  # Waveguide height
c = 5  # depth
m, n = 1, 1  # TM_mn mode (change for other modes)
epsilon = 8.854e-12  # Permittivity
mu = 4e-7 * np.pi  # Permeability

# Calculate angular frequency and wave number for default parameters
lambda_ = (3 * (10 ** 8)) / default_f
omega = 2 * np.pi * default_f
k_default = omega * np.sqrt(mu * epsilon)

# Dropdown menu for selecting m and n values
m_values = list(range(6))
n_values = list(range(6))

selected_m = tk.StringVar()
m_dropdown = ttk.Combobox(window, textvariable=selected_m, values=m_values)
m_dropdown.place(x=255, y=220, width=30)
label_m = tk.Label(window, text="m=", bg="#D9D9D9")
label_m.place(x=229, y=220)

selected_n = tk.StringVar()
n_dropdown = ttk.Combobox(window, textvariable=selected_n, values=n_values)
n_dropdown.place(x=311, y=220, width=30)
label_n = tk.Label(window, text="n=", bg="#D9D9D9")
label_n.place(x=290, y=220)

# Define a function to update the equations when m and n change
def update_equations():
    global e_x, e_y, m_x, m_y, ez_values
    e_x = ex(U, V, W)
    e_y = ey(U, V, W)
    m_x = hx(U, V, W)
    m_y = hy(U, V, W)
    ez_values = Ez(X, Y, Z, T)
    update_plot(plot_dropdown.get())

# Function to update m and n values
def update_m_n_values():
    global m, n
    m = int(selected_m.get())
    n = int(selected_n.get())

# Set default values for m and n
m_dropdown.set(str(m))
n_dropdown.set(str(n))

H = (np.pi * m / a) ** 2 + (np.pi * n / b) ** 2

def char_eqn_scaled(beta):
    return k_default ** 2 - (np.pi * m / a) ** 2 - (np.pi * n / b) ** 2 - beta ** 2

result = root(char_eqn_scaled, x0=0.1, method='hybr')
beta = result.x[0]

# Define field expressions for Ez
def Ez(x, y, z, t):
    return np.sin(m * np.pi * x / a) * np.sin(n * np.pi * y / b) * np.exp(1j * omega * t - 1 * beta * z)
A=1
def ex(x, y,z):
    return -A*(m * np.pi  / a)* np.cos(m * np.pi * x / a) * np.sin(n * np.pi * y / b) * np.exp(-1j * beta * z)
def ey(x, y,z):
    return -A*(n * np.pi / b)* np.sin(m * np.pi * x / a) * np.cos(n * np.pi * y / b) * np.exp(-1j * beta * z)
def hx(x, y,z):
    return A*(n * np.pi / b)* np.sin(m * np.pi * x / a) * np.cos(n * np.pi * y / b) * np.exp(-1j * beta * z)
def hy(x, y,z):
    return -A*(m * np.pi  / a)* np.cos(m * np.pi * x / a) * np.sin(n * np.pi * y / b) * np.exp(-1j * beta * z)

window.update_equations = update_equations

label_T = tk.Label(window, text="Time Axis(t)", bg="#D9D9D9")
label_T.place(x=19, y=110)

entry_T = ttk.Scale(window, from_=1, to=9, orient="horizontal", length=320)
entry_T.place(x=19, y=130)
entry_T.set(default_T)

label_Z = tk.Label(window, text="Space Axis(z)", bg="#D9D9D9")
label_Z.place(x=19, y=160)

entry_Z = ttk.Scale(window, from_=1, to=14, orient="horizontal", length=320)
entry_Z.place(x=19, y=180)
entry_Z.set(default_Z)

# Dropdown menu for selecting plots
plot_options = ['Hz', 'Ez', 'Ez_3d', 'contour_Ez']
selected_plot = tk.StringVar()
selected_plot.set('contour_Ez')
plot_dropdown = ttk.Combobox(window, textvariable=selected_plot, values=plot_options)
label_drop = tk.Label(window, text="Subplot=", bg="#D9D9D9")
label_drop.place(x=19, y=220)
plot_dropdown.place(x=79, y=220)

u = np.linspace(0, a, 50)
v = np.linspace(0, b, 50)
w = np.linspace(0, c, 50)
U,V,W = np.meshgrid(u, v,w)

# Create a meshgrid for the waveguide
x = np.linspace(0, a, 50)
y = np.linspace(0, b, 50)
z = np.linspace(1, 10, 50)  # Depth from 0 to 1, you can adjust this based on your requirements
t = np.linspace(1, 50, 10)
X, Y, Z, T = np.meshgrid(x, y, z, t)

def update_plot(selected_plot):

    plt.close('all')  # Close all previously opened figures

    T_new = int(entry_T.get())
    Z_new = int(entry_Z.get())

    fig=plt.figure(figsize=(9.1175, 5.3708333333))
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas_widget = canvas.get_tk_widget().place(x=379, y=14)

    if selected_plot == 'Hz':
        plt.quiver(U[:, :,0], V[:, :,0], np.real(m_x[:, :,0]), np.real(m_y[:, :,0]), scale=100, color='red')
        plt.title('Magnetic field(Hz)')
        plt.xlabel('breadth(a)')
        plt.ylabel('length(b)')
        plt.close()

    elif selected_plot == 'Ez':
        plt.quiver(U[:, :,0], V[:, :,0], np.real(e_x[:, :,0]), np.real(e_y[:, :,0]), scale=100, color='blue')
        plt.title('Electric field(Ez)')
        plt.xlabel('breadth(a)')
        plt.ylabel('length(b)')
        plt.close()

    elif selected_plot == 'Ez_3d':
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X[:, :, Z_new, T_new], Y[:, :, Z_new, T_new], np.imag(ez_values[:, :, Z_new, T_new]), cmap=cm.viridis, rstride=1, cstride=1, linewidth=0)
        ax.set_title('Electric field(3D) along z axis')
        plt.xlabel('breadth(a)')
        plt.ylabel('length(b)')
        ax.set_zlabel('Z Axia\n(depth)')
        plt.close()

    elif selected_plot == 'contour_Ez':
        plt.contourf(X[:, :, Z_new, T_new], Y[:, :, Z_new, T_new], np.real(ez_values[:, :, Z_new, T_new]), cmap=cm.cividis)
        plt.title('Magnitude plot of Electric field from end view')
        plt.xlabel('breadth(a)')
        plt.ylabel('length(b)')
        plt.close()

    fig.canvas.draw()

def plot_with_loading():

    loading_label = tk.Label(window, text="Loading...", bg="#F1F3F4",font=(14))
    loading_label.place(x=380, y=540)
    window.update()
    threading.Thread(target=update_equations()).start()
    loading_label.destroy()

# Bind the update function to dropdown selection
m_dropdown.bind("<<ComboboxSelected>>", lambda event=None: update_m_n_values())
n_dropdown.bind("<<ComboboxSelected>>", lambda event=None: update_m_n_values())
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
button_to_TE = tk.Button(window, text="TE Rectangular", command=run_program_waveguide_TE)
button_to_TE.place(x=19, y=58, width=100,height=50)

def run_program_waveguide_TMcy():
    window.destroy()
    subprocess.call(["python", "build\guiTMcy.py"])

# Button to switch to window 2
button_to_TMcy = tk.Button(window, text="TM Circular", command=run_program_waveguide_TMcy)
button_to_TMcy.place(x=122, y=58, width=100,height=50)

def run_program_waveguide_TEcy():
    window.destroy()
    subprocess.call(["python", "build\guiTEcy.py"])

button_to_TEcy = tk.Button(window, text="TE Circular", command=run_program_waveguide_TEcy)
button_to_TEcy.place(x=225, y=58, width=100,height=50)

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