import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# --- 1. Trova tutti i file .csv nella cartella corrente ---
lista_file_csv = glob.glob('amplicsv/ampli*.csv')

# Controlla se sono stati trovati file
if not lista_file_csv:
    print("Nessun file .csv trovato nella cartella.")
else:
    print(f"Trovati {len(lista_file_csv)} file .csv. Inizio l'elaborazione...")

# --- 2. Esegui un ciclo per ogni file trovato ---
for nome_file_csv in lista_file_csv:
    print(f"\n--- Elaborazione del file: {nome_file_csv} ---")
    
    # Crea il nome del file di output partendo da quello di input
    # es: 'mio_file.csv' -> 'mio_file.pdf'
    nome_base = os.path.splitext(nome_file_csv)[0]
    nome_file_grafico = f"{nome_base}.pdf"

    try:
        # Caricamento dei dati (saltando le prime 16 righe)
        dati = pd.read_csv(nome_file_csv, skiprows=16)
        
        # Se il file è vuoto o contiene solo l'header, saltalo
        if dati.empty:
            print(f"Attenzione: il file '{nome_file_csv}' è vuoto o non contiene dati. Salto al prossimo.")
            continue

        # Conversione delle unità di misura
        dati['TIME'] = dati['TIME'] * 1e9  # da s a ns
        dati['CH1'] = dati['CH1'] * 1000 # da V a mV
        dati['CH2'] = dati['CH2'] * 1000 # da V a mV
        dati['CH3'] = dati['CH3'] * 1000 # da V a mV

        # Creazione del grafico
        plt.figure(figsize=(12, 7))
        print(f"Range TIME: min={dati['TIME'].min()} ns, max={dati['TIME'].max()} ns")
        plt.xlim(-25.0, 25.0) # limit x axis from -100.0 ns to 100.0 ns
        plt.ylim(-100,100.0)
        print(f"xlim impostato a: {plt.xlim()}")        
        plt.plot(dati['TIME'], dati['CH1'], label='RPC strip')
        plt.plot(dati['TIME'], dati['CH2'], label='Scintillator trigger 1')
        plt.plot(dati['TIME'], dati['CH3'], label='Scintillator trigger 2')

        # Aggiunta di dettagli al grafico
        plt.title(f'BB5 BIS data - {nome_base}') # Titolo dinamico con il nome del file
        plt.xlabel('Time (ns)')
        plt.ylabel('Signal (mV)')
        plt.legend()
        plt.grid(True)

        # Salvataggio del grafico
        plt.savefig(nome_file_grafico)
        plt.close() # Chiude la figura per liberare memoria e non sovrapporre i grafici

        print(f"Grafico salvato con successo come '{nome_file_grafico}'")

    except Exception as e:
        print(f"Si è verificato un errore durante l'elaborazione di {nome_file_csv}: {e}")
        print("Salto al file successivo.")

print("\n--- Elaborazione completata! ---")
