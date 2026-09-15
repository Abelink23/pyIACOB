from spec import *

import pandas as pd
import tkinter as tk
from tkinter import filedialog, ttk, scrolledtext
import joblib

# print instructions message:
msg.info('The results of the selected models will be displayed in a few seconds.')
msg.info('The *** symbol highlights the class with the highest certainty obtained by the model.')
msg.info('This project is the result of the paper "Advanced ensemble techniques for the spectral classification of massive OB-type stars". For any inquiries, please contact jgonzaleze@unah.edu.hn')

# Base path relative to this script's directory
models_dir = models_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "ML-models"))

def load_models():
    global model1, model2, model3
    model1 = joblib.load(os.path.join(models_dir, "M1.joblib"))
    model2 = joblib.load(os.path.join(models_dir, "M2.joblib"))
    model3 = joblib.load(os.path.join(models_dir, "M3.joblib"))

def select_file():
    load_models()
    file_path = filedialog.askopenfilename(
        filetypes=[
            ("CSV files", "*.csv"),
            ("FITS files", "*.fits"),
            ("ASCII files", "*.ascii"),
            ("Text files", "*.txt"),
        ]
    )

    if not file_path:
        return
    if file_path:
        filename = file_path.split('/')[-1]

        if file_path.endswith("fits"):
            spectrum = spec(filename)
        else:
            spectrum = spec(filename, orig='ascii')

        spectrum.waveflux(lwl=3780.00, rwl=6849.50)
        spectrum.convolution(resol=4000)
        dlam = spectrum.dlam * len(spectrum.wave)/12279
        print(len(spectrum.wave), spectrum.wave[0], spectrum.wave[-1], dlam)
        spectrum.resamp(dlam=dlam, lwl=3780.00, rwl=6849.50, force_edges=True)
        print(len(spectrum.wave), spectrum.wave[0], spectrum.wave[-1], dlam)

        text_area.insert(tk.END, [spectrum.wave, spectrum.flux])

        models = []
        if xgb_var.get() == 1:
            models.append(model3)
        if rf_var.get() == 1:
            models.append(model1)
        if Ex_var.get() == 1:
            models.append(model2)

        predictions = []
        if len(spectrum.flux)<12279:
            dif=12279-len(spectrum.flux)
            result = [1.0 for _ in range(dif)]
            flux = np.concatenate((spectrum.flux, result))
        else:
            flux = spectrum.flux

        names = []
        for m in models:
            predictions.append(m.predict_proba(pd.DataFrame(flux).T))
            names.append(type(m).__name__)

        boxes =['O2-08','O2-08','O8.5-09.7','O8.5-09.7','B0-B2','B0-B2','B2.5-B6','B2.5-B6','B6.5-B9','B6.5-B9']
        Lumi=['III-IV-V','I-II','III-IV-V','I-II','III-IV-V','I-II','III-IV-V','I-II','III-IV-V','I-II']
        text_area.delete('1.0', tk.END)
        text_area.config(state=tk.NORMAL)
        text_area.delete('1.0', tk.END)

        text_area.insert(tk.END, f"{filename}\n")
        for name, pred in zip(names, predictions):
            text_area.insert(tk.END, f"{name}\n")
            M=np.argmax(pred)
            i=0
            for box, lu, prob in zip(boxes, Lumi, pred[0]):
                Pp=str(np.round(prob,3))
                text_area.insert(tk.END, f"\t\t{box}\t\t{lu} \t\t:\t {Pp}")
                if i==M:text_area.insert(tk.END, "   ***")
                text_area.insert(tk.END, "\n")
                i=i+1
            text_area.insert(tk.END, "\n\n")

        text_area.config(state=tk.DISABLED)


def load_image_and_text(frame):
    # Crear un label para el texto
    text_label = ttk.Label(frame, text="Texto para citar el artículo")
    text_label.pack(pady=10)

def run():
    # Crear la ventana principal
    root.title("Tool for Astronomy UNAH/IAC/UPV")
    root.resizable(False, False)

    # Crear las pestañas
    tab_control = ttk.Notebook(root)
    tab_menu = ttk.Frame(tab_control)
    tab_help = ttk.Frame(tab_control)
    tab_control.add(tab_menu, text="Main")
    tab_control.add(tab_help, text="Credits")
    tab_control.pack(expand=1, fill="both")

    # Crear el marco para los checkboxes
    frame_models = ttk.LabelFrame(tab_menu, text="Models")
    frame_models.pack(padx=10, pady=10, fill="x")

    global xgb_var, rf_var, cb_var, Ex_var
    xgb_var = tk.IntVar()
    rf_var = tk.IntVar()
    cb_var = tk.IntVar()
    Ex_var = tk.IntVar()

    ttk.Checkbutton(frame_models, text="XGB", variable=xgb_var).pack(side=tk.LEFT, padx=10)
    ttk.Checkbutton(frame_models, text="Random Forest", variable=rf_var).pack(side=tk.LEFT, padx=10)
    ttk.Checkbutton(frame_models, text="ExtraTree", variable=Ex_var).pack(side=tk.LEFT, padx=10)

    btn_select_file = ttk.Button(tab_menu, text="Select file", command=select_file)
    btn_select_file.pack(pady=10)

    global text_area
    text_area = scrolledtext.ScrolledText(tab_menu, wrap=tk.WORD, width=80, height=20)
    text_area.pack(padx=10, pady=10, fill="both", expand=True)

    text_area.config(state=tk.DISABLED)

    # Añadir la imagen y el texto en la pestaña "CREDITS"
    load_image_and_text(tab_help)

    # Ejecutar el bucle principal de la aplicación
    root.mainloop()

root = tk.Tk()
run()
