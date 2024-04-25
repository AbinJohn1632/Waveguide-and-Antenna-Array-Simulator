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
window.title("Isotropic Antennas")
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

label_Head = tk.Label(window, text="Isotropic Antennas With Feed Coefficient", bg="#D9D9D9",font=(10), relief="sunken")
label_Head.place(x=0, y=20,width=370)

# Set default values
default_de = 0
# default_N = str1,-1j])
lam_val = 1.0 

# Electric Phase (de) Slider
de_label = tk.Label(window, text="Electric Phase (de):",bg="#D9D9D9")
de_label.place(x=19, y=110)

# Fraction of Lambda Label
fraction_label = tk.Label(window, text='''Antenna Spacing (d): {:.2f} x lambda'''.format(default_de/1.0),bg="#D9D9D9")
fraction_label.place(x=19, y=165)

# Number of Antennas (N) Slider
N_label = tk.Label(window, text="Feed Coefficient",bg="#D9D9D9")
N_label.place(x=19, y=220)

de_slider = ttk.Scale(window, from_=0, to=0.9, orient="horizontal", length=320)
de_slider.set(default_de)
de_slider.place(x=19, y=135)

d_slider = ttk.Scale(window, from_=0.1, to=1.0, orient="horizontal", length=320)
d_slider.set(default_de)
d_slider.place(x=19, y=190)

N_entry = ttk.Entry(window, width=52)
# N_entry.insert(0, str(default_N))  # Set the default value
N_entry.place(x=19, y=245)

N_example = tk.Label(window, text="Example:1,1j,-1j,-1",bg="#D9D9D9")
N_example.place(x=19, y=270)


def string_to_complex_array(string_array):
    # Split the string into individual complex numbers
    complex_strings = string_array.split(',')

    # Convert each complex string to a complex number and store in a list
    complex_numbers = [complex(i) for i in complex_strings]

    # Convert the list of complex numbers to a NumPy array
    complex_array = np.array(complex_numbers, dtype=complex)
    # print(complex_array)
    return complex_array

def Afactor(w,psi):
    j = np.arange(w.size)
    # print(w)
    A = np.sum(w[j] * np.exp(j * 1j * psi[:, None]), axis=1)
    g = np.abs(A)
    return g

def para(d,de):
    """Return the power as a function of azimuthal angle, phi."""
    phi = np.linspace(0, 2*np.pi, 1000)
    psi = (2*np.pi * d / lam_val) * np.cos(phi)+de
    return psi

def get_directive_gain(g, minDdBi=-20):
    """Return the "directive gain" of the antenna array producing gain g."""
    DdBi = 10 * np.log10(g / np.max(g))
    return np.clip(DdBi, minDdBi, None)

def normal_rad_pattern(psi,N):
    pattern=(np.sin(N*psi/2)/np.sin(psi/2))/N
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
    N_val =N_entry.get()#np.array(N_entry.get())
    N_val = string_to_complex_array(N_val)  # Get the string from the entry field
    # print(N_val)

    psi = para(d_val, de_val)
    phi = np.linspace(0, 2*np.pi, 1000)
    Af=Afactor(N_val,psi)
    DdBi = get_directive_gain(Af)
    pattern = normal_rad_pattern(psi,N_val.size)

    ax.clear()
    ax.plot(phi, DdBi, linewidth=1)
    ax.set_title("Isotropic array\nAzimuthal gain Pattern\n{}".format(N_val))

    ax2.clear()
    ax2.plot(phi, pattern, linewidth=1)
    ax2.set_xlabel("phi")
    ax2.set_ylabel("pattern")
    ax2.set_title("Array factor for {} Antennas".format(N_val.size))

    # Update the fraction label
    fraction_label.config(text="Fraction of Lambda: {:.2f}".format(d_val/lam_val))

    # Update the electric phase label
    de_label.config(text="Electric Phase: {:.2f}".format(de_val))

    canvas.draw()

button_update_N = tk.Button(window, text="Plot", command=lambda value=None:update_plot())
button_update_N.place(x=125, y=300, width=100,height=50)
de_slider['command']= lambda value=None: update_plot()
d_slider['command']= lambda value=None: update_plot()

def to_home():
    window.destroy()
    subprocess.call(["python", "build\HOME.py"])

# Button to switch to window 2
button_To_home = tk.Button(window, text="HOME", command=to_home)
button_To_home.place(x=19, y=700, width=100,height=50)

def on_closing():
    window.destroy()  # Destroy the Tkinter window when closed
    exit()  # Exit the Python script, closing the terminal

def Uniform_Antennas():
    window.destroy()
    subprocess.call(["python", "build\guiAntenna1.py"])

button_to_Uniform_Antennas = tk.Button(window, text="Uniform\nAntennas", command=Uniform_Antennas)
button_to_Uniform_Antennas.place(x=122, y=700, width=100,height=50)

# print(psi,phi)
window.protocol("WM_DELETE_WINDOW", on_closing)
window.resizable(0, 0)
window.mainloop()