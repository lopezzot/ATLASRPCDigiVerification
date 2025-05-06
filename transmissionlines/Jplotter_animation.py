import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import re # Import regular expressions for parsing filenames

# --- Configurazione ---
# Modifica questi valori secondo le tue necessità
FILE_PATTERN = "J*.txt" # Pattern per trovare i file (es. rpc_signal_t*ns.txt)
                                   # Assicurati che il prefisso corrisponda!
OUTPUT_FILENAME = "J_animation.gif" # Nome file output (.gif o .mp4)
INTERVAL_MS = 100  # Millisecondi tra un frame e l'altro (velocità animazione)

# Limiti asse Y: None per auto-determinazione, oppure una tupla (ymin, ymax)
# Es: Y_LIMITS = (-0.01, 0.01) # Imposta manualmente se necessario
Y_LIMITS = None

# --- Funzione per Estrarre il Tempo dal Nome File ---
def extract_time_from_filename(filename):
    """Estrae il tempo in ns da nomi file tipo 'prefix_t12.345ns.txt'."""
    # Cerca un pattern che corrisponde a '_t' seguito da un numero (anche negativo, con decimali)
    # e poi 'ns.txt' alla fine del nome base del file.
    match = re.search(r'_t(-?\d+\.?\d*)ns\.txt$', os.path.basename(filename))
    if match:
        try:
            # Estrae il numero catturato (gruppo 1) e lo converte in float
            return float(match.group(1))
        except ValueError:
            # Se la conversione fallisce, restituisce None
            return None
    # Se il pattern non matcha, restituisce None
    return None

# --- Trova e Ordina i File ---
print(f"Ricerca file corrispondenti a: '{FILE_PATTERN}'")
file_list = glob.glob(FILE_PATTERN)

if not file_list:
    print(f"Errore: Nessun file trovato per il pattern '{FILE_PATTERN}'. Controlla il percorso e il pattern.")
    exit()

# Crea una lista di dizionari per associare file e tempo, gestendo errori
file_data = []
for f in file_list:
    time_ns = extract_time_from_filename(f)
    if time_ns is not None:
        file_data.append({'filename': f, 'time_ns': time_ns})
    else:
        print(f"Attenzione: Impossibile estrarre il tempo da '{f}'. File saltato.")

# Controlla se abbiamo file validi dopo il parsing
if not file_data:
    print("Errore: Nessun file valido trovato dopo il parsing dei nomi. Controlla il formato dei nomi file.")
    exit()

# Ordina la lista basandosi sul tempo estratto
file_data.sort(key=lambda x: x['time_ns'])

# Liste ordinate di nomi file e tempi
sorted_files = [item['filename'] for item in file_data]
times_ns = [item['time_ns'] for item in file_data]
num_frames = len(sorted_files)

print(f"Trovati {num_frames} file validi da animare.")

# --- Funzione per Leggere i Dati da un File ---
def read_data(filename):
    """Legge dati x, V da un file txt, saltando header con #."""
    try:
        # Usa comments='#' per ignorare le righe di header
        data = np.loadtxt(filename, comments='#')
        # Assicura che ci siano almeno 2 colonne
        if data.ndim == 2 and data.shape[1] >= 2:
            return data[:, 0], data[:, 1] # Colonna 0: x, Colonna 1: V
        elif data.ndim == 1 and data.shape[0] >= 2: # Caso possibile se c'è solo un punto dati?
             print(f"Attenzione: Formato dati 1D inatteso in {filename}. Provo a interpretare.")
             # Potrebbe essere necessario adattare se il formato è diverso
             return None, None # Meglio saltare se incerto
        else:
            print(f"Attenzione: Formato dati non valido (shape: {data.shape}) in {filename}. File saltato.")
            return None, None
    except Exception as e:
        print(f"Errore durante la lettura di {filename}: {e}")
        return None, None

# --- Determina Limiti Assi (Consigliato per stabilità visualizzazione) ---
x_min, x_max = 0, 0
y_min_global, y_max_global = np.inf, -np.inf

# Leggi il primo file per ottenere la struttura dell'asse x
x_data_first, v_data_first = read_data(sorted_files[0])
if x_data_first is None:
    print("Errore: Impossibile leggere i dati dal primo file per inizializzare il plot.")
    exit()
x_min, x_max = np.min(x_data_first), np.max(x_data_first)
y_min_global = min(y_min_global, np.min(v_data_first))
y_max_global = max(y_max_global, np.max(v_data_first))

# Se Y_LIMITS non sono impostati manualmente, leggiamo gli altri file
if Y_LIMITS is None:
    print("Determinazione automatica limiti asse Y (potrebbe richiedere tempo)...")
    for filename in sorted_files[1:]: # Inizia dal secondo file
        _, v_data = read_data(filename)
        if v_data is not None:
            y_min_global = min(y_min_global, np.min(v_data))
            y_max_global = max(y_max_global, np.max(v_data))

    # Aggiungi un po' di padding ai limiti y auto-determinati
    if np.isinf(y_min_global) or np.isinf(y_max_global):
        print("Attenzione: Impossibile determinare limiti Y validi. Uso default (-1, 1).")
        final_y_limits = (-1, 1)
    else:
        padding = (y_max_global - y_min_global) * 0.1
        padding = max(padding, 0.01) # Assicura un minimo padding se i valori sono quasi costanti
        final_y_limits = (y_min_global - padding, y_max_global + padding)
    print(f"Limiti Y calcolati: ({final_y_limits[0]:.3f}, {final_y_limits[1]:.3f})")
else:
    final_y_limits = Y_LIMITS # Usa i limiti impostati manualmente
    print(f"Uso limiti Y manuali: {final_y_limits}")


# --- Setup del Plot ---
fig, ax = plt.subplots()
# Inizializza una linea vuota che verrà aggiornata
# usiamo i dati x del primo file come riferimento per l'asse x
line, = ax.plot(x_data_first, np.zeros_like(x_data_first), lw=2) # Inizia con V=0
ax.set_xlabel("Strip position (m)")
ax.set_ylabel("Current density (A/m)")
ax.set_title("RPC Current density from a single spike")
ax.set_xlim(x_min, x_max)
ax.set_ylim(final_y_limits)
ax.grid(True)
# Aggiungi testo per mostrare il tempo corrente sull'animazione
time_text = ax.text(0.05, 0.90, '', transform=ax.transAxes, fontsize=10,
                    bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))

# --- Funzione di Inizializzazione (per blitting) ---
# Reimposta la linea e il testo all'inizio o quando la figura viene ridisegnata
def init():
    line.set_ydata(np.zeros_like(x_data_first)) # Reimposta a zero
    time_text.set_text('')
    return line, time_text

# --- Funzione di Aggiornamento (chiamata per ogni frame) ---
def update(frame_index):
    """Aggiorna i dati della linea per il frame corrente."""
    filename = sorted_files[frame_index]
    time_ns = times_ns[frame_index]

    x, v = read_data(filename)

    if x is not None and v is not None:
        # Assumiamo che l'asse x non cambi tra i file
        # Se potesse cambiare, dovremmo usare line.set_data(x, v)
        line.set_ydata(v)
        time_text.set_text(f'Time = {time_ns:.3f} ns')
    else:
        # Se c'è un errore di lettura, mostra linea vuota e segnalalo
        line.set_ydata(np.zeros_like(x_data_first)) # Linea a zero
        time_text.set_text(f'Time = {time_ns:.3f} ns\n(Errore Lettura File!)')

    # Ritorna una tupla degli oggetti grafici modificati (necessario per blit=True)
    return line, time_text

# --- Crea l'Oggetto Animazione ---
# fig: la figura matplotlib
# update: la funzione chiamata ad ogni frame
# frames: il numero totale di frame (numero di file)
# init_func: la funzione per inizializzare il plot
# blit=True: ottimizzazione che ridisegna solo le parti cambiate (richiede init_func e che update ritorni gli artisti modificati)
# interval: intervallo tra frame in ms
# repeat=False: non ripetere l'animazione alla fine
print("Creazione animazione...")
ani = animation.FuncAnimation(fig, update, frames=num_frames,
                            init_func=init, blit=True, interval=INTERVAL_MS, repeat=False)

# --- Salva o Mostra l'Animazione ---
try:
    # Determina lo 'writer' in base all'estensione del file e a cosa è installato
    if OUTPUT_FILENAME.lower().endswith(".gif"):
        print(f"Salvataggio animazione come GIF in '{OUTPUT_FILENAME}' (potrebbe richiedere Pillow)...")
        # Pillow è spesso una dipendenza di matplotlib o facile da installare
        ani.save(OUTPUT_FILENAME, writer='pillow', fps=1000 / INTERVAL_MS)
    elif OUTPUT_FILENAME.lower().endswith(".mp4"):
        print(f"Salvataggio animazione come MP4 in '{OUTPUT_FILENAME}' (richiede FFmpeg)...")
        # Richiede che ffmpeg sia installato e nel PATH di sistema
        ani.save(OUTPUT_FILENAME, writer='ffmpeg', fps=1000 / INTERVAL_MS)
    else:
        # Formato non riconosciuto, prova a salvare come GIF
        print(f"Estensione file non riconosciuta. Provo a salvare come GIF in '{OUTPUT_FILENAME}.gif'...")
        ani.save(OUTPUT_FILENAME + ".gif", writer='pillow', fps=1000 / INTERVAL_MS)
    print("Salvataggio completato.")

except Exception as e:
    print("\n---------------------------------------------------------")
    print(f"Errore durante il salvataggio dell'animazione: {e}")
    print("Possibili cause:")
    print(" - Lo 'writer' necessario non è installato (es. Pillow per GIF, FFmpeg per MP4).")
    print("   Prova a installare: 'pip install Pillow' o installa FFmpeg dal suo sito.")
    print(" - Problemi di permessi di scrittura nella cartella.")
    print("Mostro l'animazione direttamente (chiudi la finestra per terminare).")
    print("---------------------------------------------------------")
    plt.show()

# In alternativa, per mostrare direttamente l'animazione senza salvarla:
# print("Mostro l'animazione (chiudi la finestra per terminare)...")
# plt.show()
