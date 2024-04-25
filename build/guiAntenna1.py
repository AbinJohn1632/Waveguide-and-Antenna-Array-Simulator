from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage, filedialog,ttk, Label, Scrollbar
import threading
import subprocess

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"E:\python\machine_learning\college_project\Waveguide_Simulator_sem6\build\assets\frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

window = tk.Tk()
window.title("Uniform linear array Antenna")
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

label_Head = tk.Label(window, text="Uniform Linear Array Antenna", bg="#D9D9D9",font=(10), relief="sunken")
label_Head.place(x=0, y=20,width=370)

# Set default values
default_de = 0
default_N = 1
lam_val = 1.0  # Set lambda=1

# Electric Phase (de) Slider
de_label = tk.Label(window, text="Electric Phase (de):",bg="#D9D9D9")
de_label.place(x=19, y=110)

# Fraction of Lambda Label
fraction_label = tk.Label(window, text='''Antenna Spacing (d): {:.2f} x lambda'''.format(default_de/1.0),bg="#D9D9D9")
fraction_label.place(x=19, y=160)

# Number of Antennas (N) Slider
N_label = tk.Label(window, text="Number of Antennas (N):",bg="#D9D9D9")
N_label.place(x=19, y=210)

de_slider = ttk.Scale(window, from_=0, to=0.9, orient="horizontal", length=320)
de_slider.set(default_de)
de_slider.place(x=19, y=130)

d_slider = ttk.Scale(window, from_=0.1, to=1.0, orient="horizontal", length=320)
d_slider.set(default_de)
d_slider.place(x=19, y=180)

N_slider = ttk.Scale(window, from_=1, to=20, orient="horizontal", length=320)
N_slider.set(default_N)
N_slider.place(x=19, y=230)

def para(d, w, de):
    lam_val = 1.0
    phi = np.linspace(0, 2*np.pi, 1000)
    psi = (2*np.pi * d / lam_val) * np.cos(phi) + de
    return psi, phi

def normal_rad_pattern(psi, N):
    pattern = (np.sin(N*psi/2)/np.sin(psi/2))/N
    return pattern

def update_plot(event=None):
    plt.close('all')

    fig=plt.figure(figsize=(9.1175, 5.3708333333))
    ax = fig.add_subplot(121, polar=True)
    ax2 = fig.add_subplot(122)

    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget = canvas.get_tk_widget().place(x=379, y=14)

    de_val = float(de_slider.get())
    d_val = float(d_slider.get())
    N_val = int(N_slider.get())

    psi, phi = para(d_val, N_val, de_val)
    pattern = normal_rad_pattern(psi, N_val)

    ax.clear()
    ax.plot(phi, pattern, linewidth=1)
    ax.set_title("Uniform linear array\nAzimuthal gain Pattern N={}".format(N_val))

    ax2.clear()
    ax2.plot(phi, pattern, linewidth=1)
    ax2.set_xlabel("phi")
    ax2.set_ylabel("pattern")
    ax2.set_title("Uniform linear array\nArray factor for N={}".format(N_val))

    # Update the fraction label
    fraction_label.config(text="Antenna Spacing (d): {:.2f} x lambda".format(d_val/lam_val))

    # Update the electric phase label
    de_label.config(text="Electric Phase: {:.2f}".format(de_val))

    canvas.draw()

de_slider['command']= lambda value=None: update_plot()
d_slider['command']= lambda value=None: update_plot()
N_slider['command']= lambda value=None: update_plot()

button = tk.Button(window, text="Plot", command=lambda value=None:update_plot())
button.place(x=125, y=270, width=100,height=50)

def to_home():
    window.destroy()
    subprocess.call(["python", "build\HOME.py"])

# Button to switch to window 2
button_To_home = tk.Button(window, text="HOME", command=to_home)
button_To_home.place(x=19, y=700, width=100,height=50)

def Isotropic_Antennas():
    window.destroy()
    subprocess.call(["python", "build\guiAntenna2.py"])

button_to_Isotropic_Antennas = tk.Button(window, text="Isotropic\nAntennas", command=Isotropic_Antennas)
button_to_Isotropic_Antennas.place(x=122, y=700, width=100,height=50)

def on_closing():
    window.destroy()  # Destroy the Tkinter window when closed
    exit()  # Exit the Python script, closing the terminal

window.protocol("WM_DELETE_WINDOW", on_closing)
# Run the Tkinter event loop
window.resizable(0, 0)
window.mainloop()