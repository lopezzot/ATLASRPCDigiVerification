#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

// Function to save profile of a vector
void save_profile(const std::vector<double>& V, double time, int N, double dx, const std::string& prefix = "output") {
    std::ostringstream filename;
    filename << prefix << "_t" << std::fixed << std::setprecision(5) << time * 1e9 << "ns.txt";
    std::ofstream out(filename.str());
     if (!out) {
        std::cerr << "Error while opening the file: " << filename.str() << std::endl;
        return;
    }
    out << std::fixed << std::setprecision(6);
    // Save position and value
    out << "# x (m)\tValue\n";
    for (int i = 0; i < N; ++i) {
        out << i * dx << "\t" << V[i] << "\n";
    }
}
void save_timefactor(const std::vector<double>& V, double dt, const std::string& prefix = "timefactor") {
    std::ostringstream filename;
    filename << prefix << ".txt";
    std::ofstream out(filename.str());
     if (!out) {
        std::cerr << "Error while opening the file: " << filename.str() << std::endl;
        return;
    }
    out << std::fixed << std::setprecision(6);
    // Save position and value
    out << "# t (ns)\tValue\n";
    for (int i = 0; i < V.size(); ++i) {
        out << (i+0.5)*dt*1e9 << "\t" << V[i] << "\n";
    }
}

#define wLOSS

int main() {
    // --- Parameters for grid and transmission line ---
    const int N = 500;        // Numero di punti griglia
    const double length = 2.0; // Lunghezza fisica striscia (m)
    const double dx = length / (N - 1); // Passo spaziale (m) (segmento di linea)
    const double L = 2.5e-7;   // Induttanza per unità di lunghezza (H/m)
    const double C = 1e-10;  // Capacità per unità di lunghezza (F/m)
    #ifdef wLOSS
    const double R = 1.0;    // Resistenza per unità di lunghezza (Ohm/m) - VALORE ESEMPIO!
    #endif
    const double Z0 = std::sqrt(L / C); // Impedenza caratteristica (Ohm)
    const double R_L = Z0; // <<< NUOVA RIGA: Resistenza di carico ai bordi [Ohm] - Assumiamo adattata
                         //
    // --- Parametri Temporali e Stabilità ---
    const double v = 1.0 / std::sqrt(L * C); // Velocità di propagazione (m/s)
    const double dt = 0.95 * dx / v;         // Passo temporale (s) - Fattore CFL 0.95
    const double mur_const = (v * dt - dx) / (v * dt + dx); // Costante BC Mur

    // --- Parametri Simulazione ---
    const double T = 2.0 * (length / v); // Max time of computing
    const int steps = static_cast<int>(T / dt);
    const int snapshot_interval = steps / 20; // Salva circa 100 snapshot
    const double output_precision = 1e-4; // Soglia per considerare V trascurabile e terminare programma

    // --- Parametri Sorgente RPC Realistica ---
    // !!! QUESTI VALORI SONO ESEMPI - USA VALORI SPECIFICI PER ATLAS RPC !!!
    const double x_source = length / 2.0; // Posizione sorgente (es. 1/3 della lunghezza) (m)
    const double sigma_x = 20e-3;         // Larghezza spaziale sorgente (es. 3 mm) (m)
    const double t_start = 0.;//10e-9;        // Tempo inizio impulso sorgente (s) (es. 10 ns)
    const double tau = 0.5e-9;//2.5e-9;           // Costante di tempo impulso (s) (es. 2.5 ns)
    const double J_peak = -5e-3;         // Picco densità corrente [A/m] (SEGNO NEGATIVO = carica indotta) - VALORE DA CALIBRARE!

    std::cout << "--- Parametri Simulazione RPC ---" << std::endl;
    std::cout << " N = " << N << ", length = " << length << " m, dx = " << dx << " m" << std::endl;
    std::cout << " L = " << L << " H/m, C = " << C << " F/m, Z0 = " << Z0 << " Ohm" << std::endl;
    std::cout << " v = " << v << " m/s, dt = " << dt * 1e9 << " ns" << std::endl;
    std::cout << " T_max = " << T * 1e9 << " ns, steps = " << steps << std::endl;
    std::cout << " BC Mur Const = " << mur_const << std::endl;
    std::cout << "--- Parametri Sorgente ---" << std::endl;
    std::cout << " x_source = " << x_source << " m, sigma_x = " << sigma_x << " m" << std::endl;
    std::cout << " t_start = " << t_start * 1e9 << " ns, tau = " << tau * 1e9 << " ns" << std::endl;
    std::cout << " J_peak = " << J_peak << " A/m" << std::endl;
    std::cout << "----------------------------------" << std::endl;

    #ifdef wLOSS
    std::cout << "--- Including resistivity --" << std::endl;
    std::cout << " R = " << R << " Ohm/m" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    #endif

    // Vettori Tensione e Corrente - Inizializzati a zero
    std::vector<double> V(N, 0.0), V_new(N, 0.0);
    std::vector<double> I(N, 0.0), I_new(N, 0.0); // I[i] è tra i e i+1
    std::vector<double> J(N, 0.0); // input densità di corrente [A/m] dalla scarica degli rpc
    std::vector<double> time_factor(0.0);

    // Variabili per rilevamento bordi
    //const double threshold = std::abs(J_peak * Z0 * 0.1); // probably an error, (A/m)*ohm)=V/m cannot compare with V Soglia dinamica (es. 10% ampiezza stimata)
    const double threshold = 0.001; // V
    bool reached_left = false, reached_right = false;
    double time_left = -1.0, time_right = -1.0; // Inizializza a -1

    std::cout << "--- Threshold to detect signal at borders ---" << std::endl;
    std::cout << " threshold : " << threshold << " V " << std::endl;
    std::cout << "----------------------------------" << std::endl;

    // Costanti moltiplicative per aggiornamento FDTD
    const double const_I = dt / (L * dx); // dt/(L*dx)
    const double const_V = dt / (C * dx); // dt/(C*dx)
    const double const_Source = dt / C;   // dt/C per termine sorgente

    #ifdef wLOSS
    const double denom_I = 1.0 + R * dt / (2.0 * L);
    const double cI1 = (1.0 - R * dt / (2.0 * L)) / denom_I;
    const double cI2 = (dt / (L * dx)) / denom_I;

     // <<< NUOVO BLOCCO: Costanti per BC Resistiva (G=0) >>>
    double cV1_bc = 0.0, cV2_bc = 0.0;
    const double R_load = R_L; // Usa R_L definita sopra
    if (std::abs(R_load) > 1e-18) { // Evita divisione per zero se R_L=0 (cortocircuito)
        const double c1_bc = C * dx / dt;
        const double c2_bc = 1.0 / (2.0 * R_load);
        const double den_bc = c1_bc + c2_bc;
        if (std::abs(den_bc) < 1e-18) {
             std::cerr << "Errore: den_bc (BC coeff) vicino a zero." << std::endl; return 1;
        }
        cV1_bc = (c1_bc - c2_bc) / den_bc;
        cV2_bc = 1.0 / den_bc;
    } else {
        // Se R_L == 0 (cortocircuito), V al bordo è forzato a 0
        cV1_bc = 0.0;
        cV2_bc = 0.0; // V sarà 0
    }
    #endif

    // --- Loop Temporale FDTD ---
    for (int step = 0; step < steps; ++step) {

        // Aggiorna Corrente I (da i=0 a N-2) -> I_new è I^{n+1/2}
        for (int i = 0; i < N - 1; ++i) {
            #ifndef wLOSS
            I_new[i] = I[i] - const_I * (V[i + 1] - V[i]);
            #endif
            #ifdef wLOSS
            I_new[i] = cI1 * I[i] - cI2 * (V[i + 1] - V[i]);
            #endif
        }
        // Aggiorna Tensione V (punti interni: da i=1 a N-2) -> V_new è V^{n+1} parziale
        // V^{n+1}_i = V^n_i - (dt/C*dx) * (I^{n+1/2}_{i+1/2} - I^{n+1/2}_{i-1/2})
        for (int i = 1; i < N - 1; ++i) {
            V_new[i] = V[i] - const_V * (I_new[i] - I_new[i - 1]);
        }

        // Aggiungi Termine Sorgente J(x,t) a V_new (punti interni)
        double current_sim_time = (step + 0.5) * dt; // Tempo centrato per sorgente
        double time_relative = current_sim_time - t_start;

        // Calcola sorgente solo se nel periodo attivo
        if (time_relative > 0) {
            // Forma temporale: (t'/tau) * exp(-t'/tau) normalizzata al picco
            // Modifichiamo per avere un integrale definito (carica totale per unità di lunghezza)
            // Usiamo Q_density * (1/tau) * exp(-t'/tau) dove Q_density = J_peak * tau
            // Oppure usiamo J_peak * (t'/tau^2) * exp(-t'/tau)
            double time_norm = time_relative / tau;
            double temporal_factor = 0.0;
            if (time_norm > 0){
            // Scaliamo J_peak per rappresentare la densità di corrente A/m al picco temporale
            // Il picco di (t'/tau)exp(-t'/tau) è 1/e a t'=tau.
            // Il picco di (t'/tau^2)exp(-t'/tau) è 1/(e*tau) a t'=tau
            // Se J_peak è A/m, usiamo:
             temporal_factor = time_norm * std::exp(-time_norm); // Forma proporzionale alla corrente
             // Normalizza in modo che il picco della corrente temporale sia 1 a t'=tau
             // il picco di (t'/tau)*exp(-t'/tau) è 1/e
             temporal_factor *= std::exp(1.0); // Ora il picco è 1
             time_factor.push_back(temporal_factor);
            }
            for (int i = 1; i < N - 1; ++i) {
                double x = i * dx;
                double spatial_factor = std::exp(-std::pow((x - x_source), 2) / (2 * sigma_x * sigma_x));
                if(step<100 && step % 5 ==0){
                  J[i] = J_peak * spatial_factor * temporal_factor;
                  save_profile(J, (step+0.5)*dt, N, dx, "J"); // Passa N e dx
                }
                double source_current_density = J_peak * spatial_factor * temporal_factor;

                // Applica a V: V_new[i] -= (dt / C) * J(x,t)
                V_new[i] -= const_Source * source_current_density;
            }
        }

        #ifndef wLOSS
        // Applica BC assorbenti di Mur per V[0] e V[N-1] (sovrascrive V_new ai bordi)
        // Usa V (valori a n) e V_new (valori a n+1 dei punti interni già calcolati)
        V_new[0] = V[1] + mur_const * (V_new[1] - V[0]);

        if (N > 1) {
            V_new[N - 1] = V[N - 2] + mur_const * (V_new[N - 2] - V[N - 1]);
        } else {
             V_new[0] = 0; // Caso N=1
        } // fine boundary condition*/
        #endif
 
        #ifdef wLOSS
        // <<< INIZIO BLOCCO BC RESISTIVE (sostituisce il blocco Mur) >>>
        if (std::abs(R_L) < 1e-18) { // Cortocircuito
          if (N > 0) V_new[0] = 0.0;
          if (N > 1) V_new[N - 1] = 0.0;
        }
        else { // Carico resistivo finito
            // Bordo Sinistro (V_new[0]) - Usa V[0], I_new[0]
            if (N > 1) {
             V_new[0] = cV1_bc * V[0] - cV2_bc * I_new[0];
            } else if (N == 1) { V_new[0] = 0; } // Non dovrebbe succedere con N=500

          // Bordo Destro (V_new[N-1]) - Usa V[N-1], I_new[N-2]
          if (N > 1) {
          V_new[N - 1] = cV1_bc * V[N - 1] + cV2_bc * I_new[N - 2];
          }
        }
        // <<< FINE BLOCCO BC RESISTIVE >>>
        #endif

        // Tempo effettivo alla fine dello step
        double t_output = (step + 1) * dt;
        double t_output_forcurrent = (step+0.5) * dt;

        // Controlla arrivo ai bordi (usa V_new appena calcolato)
        if (!reached_left && std::abs(V_new[0]) >= threshold) { // Usa abs per soglia
            reached_left = true;
            time_left = t_output;
            // std::cout << "-> Segnale raggiunto bordo sinistro a t = " << t_output*1e9 << " ns" << std::endl;
        }

        if (!reached_right && std::abs(V_new[N - 1]) >= threshold) { // Usa abs per soglia
            reached_right = true;
            time_right = t_output;
            // std::cout << "-> Segnale raggiunto bordo destro a t = " << t_output*1e9 << " ns" << std::endl;
        }

        // Scrivi su file snapshot
        if (t_output < 1.5*tau && step % 2 == 0){ 
           save_profile(V_new, t_output, N, dx, "rpc_signal"); // Passa N e dx
           save_profile(I_new, t_output_forcurrent, N, dx, "I"); // Passa N e dx
        }
        if ((t_output > 1.5*tau && step % snapshot_interval == 0) || step == steps -1 ) { // Salva anche l'ultimo
           // std::cout << "Salvataggio snapshot a t = " << t_output*1e9 << " ns (step " << step << ")" << std::endl;
           save_profile(V_new, t_output, N, dx, "rpc_signal"); // Passa N e dx
           save_profile(I_new, t_output_forcurrent, N, dx, "I"); // Passa N e dx
        }

        // Aggiorna stato per il prossimo ciclo
        V.swap(V_new);
        I.swap(I_new);

        // Condizione di uscita anticipata (opzionale)
        if (reached_left && reached_right && step > steps / 2 ) {
             double max_V_abs = 0.0;
             for(double val : V) { max_V_abs = std::max(max_V_abs, std::abs(val)); }
             if (max_V_abs < output_precision) {
                  std::cout << "\nSegnale uscito e tensione residua bassa a t = " << t_output*1e9 << " ns. Interrompo." << std::endl;
                  // Salva ultimo stato se non già salvato
                   if (step % snapshot_interval != 0) {
                       save_profile(V, t_output, N, dx, "rpc_signal");
                       save_profile(I, t_output_forcurrent, N, dx, "I");
                   }
                  break;
             }
        }

    } // Fine loop temporale

    save_timefactor(time_factor, dt);

    std::cout << "\n--- Risultati Rilevamento Bordi ---" << std::endl;
    if (time_left > 0)
       std::cout << "Tempo di arrivo a sinistra: " << time_left * 1e9 << " ns\n";
    else
       std::cout << "Segnale non rilevato a sinistra (soglia=" << threshold <<").\n";

    if (time_right > 0)
       std::cout << "Tempo di arrivo a destra: " << time_right * 1e9 << " ns\n";
    else
       std::cout << "Segnale non rilevato a destra (soglia=" << threshold <<").\n";

    std::cout << "Velocità teorica di propagazione: " << v << " m/s\n";

    return 0;
}
