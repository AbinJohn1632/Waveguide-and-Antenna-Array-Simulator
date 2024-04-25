from pathlib import Path
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage
import tkinter as tk

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"E:\python\machine_learning\college_project\Waveguide_Simulator_sem6\build\assets\frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

window = Tk()

window.geometry("1305x770")
window.configure(bg = "#FFFFFF")

canvas = Canvas(
    window,
    bg = "#FFFFFF",
    height = 770,
    width = 1305,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge"
)

canvas.place(x = 0, y = 0)
image_image_1 = PhotoImage(
    file=relative_to_assets("image_1.png"))
image_1 = canvas.create_image(
    652.0,
    425.0,
    image=image_image_1
)

import subprocess

def run_program_waveguide():
    window.withdraw()
    subprocess.call(["python", "build\guiTErec.py"])

button_image_1 = PhotoImage(
    file=relative_to_assets("button_2.png"))

button_1 = tk.Button(window,
    image=button_image_1,
    borderwidth=0,
    highlightthickness=0,
    command=run_program_waveguide,
    relief="flat"
)

button_1.place(
    x=511.0,
    y=374,
    width=285.0,
    height=45.0
)

def run_program_antenna():
    window.withdraw()
    subprocess.call(["python", "build\guiAntenna1.py"])

button_image_2 = PhotoImage(
    file=relative_to_assets("button_1.png"))

button_2 = tk.Button(
    image=button_image_2,
    borderwidth=0,
    highlightthickness=0,
    command=run_program_antenna,
    relief="flat"
)

button_2.place(
    x=511.0,
    y=450,
    width=285.0,
    height=45.0
)

def on_closing():
    window.destroy()  # Destroy the Tkinter window when closed
    exit()  # Exit the Python script, closing the terminal

window.protocol("WM_DELETE_WINDOW", on_closing)

window.resizable(False, False)
window.mainloop()
