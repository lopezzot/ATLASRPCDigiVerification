#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <random>

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
// Funzione per salvare l'output del Charge Sensitive Amplifier (esempio)
void save_csa_output(const std::vector<double>& times,
                     const std::vector<double>& csa_left,
                     const std::vector<double>& csa_right,
                     const std::string& filename = "csa_output.txt") {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error opening CSA output file: " << filename << std::endl;
        return;
    }
    out << std::fixed << std::setprecision(6);
    out << "# Time (ns)\tV_CSA_Left (V)\tV_CSA_Right (V)\n";
    for (size_t i = 0; i < times.size(); ++i) {
        out << times[i] * 1e9 << "\t"
            << (i < csa_left.size() ? csa_left[i] : 0.0) << "\t"
            << (i < csa_right.size() ? csa_right[i] : 0.0) << "\n";
    }
    //std::cout << "CSA output saved to " << filename << std::endl;
}

struct rpcoutput{

    // metadata needed as input
    double length; // m
    int N;
    double R; // ohm/m
    double tau; //s
    double threshold; // V
    double Jpeak; // A/m
    double SignalJitter; // s
    double TDCBinSize; // s

    void add_metadata(double aLength, int aN, double aR, double aTau, double aThreshold, double aJpeak, double aSignalJitter, double aTDCBinSize){
        length=aLength;
        N=aN;
        R=aR;
        tau=aTau;
        threshold=aThreshold;
        Jpeak = aJpeak;
        SignalJitter = aSignalJitter;
        TDCBinSize = aTDCBinSize;
    };
    
    // data from simulation output
    double time_left; // s
    double time_right; // s
    double tot_left; // time-over-threshold left (s)
    double tot_right; // time-over-threshold right (s)

    void add_output(double aTL, double aTR, double aTOTleft, double aTOTright){
        time_left = aTL;
        time_right = aTR;
        tot_left = aTOTleft;
        tot_right = aTOTright;
    };

    // fields and methods for file creation
    std::string fileName = "output_rpc.txt";
    void set_filename(std::string& aString){ fileName = aString; };

    void rpcoutput_tofile(bool append_output){
        std::ostringstream filename;
        filename << fileName;
        std::ofstream out;
        if(!append_output) out.open(filename.str());
        else out.open(filename.str(), std::ios::app);
        if (!out) {
            std::cerr << "Error while opening the file: " << filename.str() << std::endl;
            return;
        }
        out << std::fixed << std::setprecision(6);

        // Save both input and output fields
        if(!append_output) out << "length (m)\t N\t R (ohm/m)\t tau (ns)\t threshold (V)\t time_left (ns)\t time_right (ns)\t tot_left (ns)\t tot_right (ns)\t Jpeak (A/m)\t Jitter (ns)\t TDCBin (s) \n";
        out << length << "\t" << N << "\t" << R << "\t" << tau*1e9 << "\t" << threshold << "\t" << time_left*1e9 << "\t" << time_right*1e9 << "\t" << tot_left*1e9 << "\t" << tot_right*1e9 << "\t" << Jpeak << "\t" << SignalJitter*1e9 << "\t" << TDCBinSize*1e9 <<"\n";
    }
};

#define wLOSS
//#define MURBC

void process_rpc_signal(double aLength, int aN, [[maybe_unused]] double aR, double aTau, double aThreshold, double aJpeak, double aSignaljitter, double aTDCbinsize, std::string output_name = "output_rpc.txt", bool append_output = false) {
    // --- Parameters for grid and transmission line ---
    const int N = aN;        // Numero di punti griglia
    const double length = aLength; // Lunghezza fisica striscia (m)
    const double dx = length / (N - 1); // Passo spaziale (m) (segmento di linea)
    //const double L = 2.5e-7;   // Induttanza per unità di lunghezza (H/m)
    const double L = 9.01e-8; // Induttanza calcolata per avere impedenza 18 ohm e velocità 20cm/ns (H/m)
    //const double C = 1e-10;  // Capacità per unità di lunghezza (F/m)
    const double C = 2.78e-10; // Capacità calcolata per avere impedenza 18 ohm e velocità 20cm/ns (F/m)
    #ifdef wLOSS
    const double R = aR;    // Resistenza per unità di lunghezza (Ohm/m) - VALORE ESEMPIO!
    #endif
    const double Z0 = std::sqrt(L / C); // Impedenza caratteristica (Ohm)
    const double R_L = Z0; // <<< NUOVA RIGA: Resistenza di carico ai bordi [Ohm] - Assumiamo adattata

    // --- Parametri Temporali e Stabilità ---
    const double v = 1.0 / std::sqrt(L * C); // Velocità di propagazione (m/s)
    const double dt = 0.95 * dx / v;         // Passo temporale (s) - Fattore CFL 0.95
    const double mur_const = (v * dt - dx) / (v * dt + dx); // Costante BC Mur

    // --- Parametri Simulazione ---
    double T = 0.;
    if(length<0.8) T = 40.0 * (length / v); // Max time of computing
    else T = 5.0 * (length/v);
    const int steps = static_cast<int>(T / dt);
    int snapshot_interval = steps / 20; // Salva circa 100 snapshot
    const double output_precision = 1e-4; // Soglia per considerare V trascurabile e terminare programma

    // --- Parametri Sorgente RPC Realistica ---
    // !!! QUESTI VALORI SONO ESEMPI - USA VALORI SPECIFICI PER ATLAS RPC !!!
    const double x_source = length / 2.0; // Posizione sorgente (es. 1/3 della lunghezza) (m)
    const double sigma_x = 20e-3;         // Larghezza spaziale sorgente (es. 3 mm) (m)
    const double t_start = 0.;//10e-9;        // Tempo inizio impulso sorgente (s) (es. 10 ns)
    const double tau = aTau;//0.5e-9;//2.5e-9;           // Costante di tempo impulso (s) (es. 2.5 ns)
    const double J_peak = aJpeak;//-7e-3;         // Picco densità corrente [A/m] (SEGNO NEGATIVO = carica indotta) - VALORE DA CALIBRARE!
    // calcolo carica iniettata sulla strip
    constexpr double e = 2.718281828459045;
    const double ChargeOnStrip = -1.*J_peak * std::sqrt(2*M_PI*sigma_x*sigma_x) * e*tau;

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
    std::cout << " Charge on strip = " << ChargeOnStrip * 1e15 << " fC " << std::endl;
    std::cout << "----------------------------------" << std::endl;

    #ifdef wLOSS
    std::cout << "--- Including resistivity --" << std::endl;
    std::cout << " R = " << R << " Ohm/m" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    #endif

    // Vettori Tensione e Corrente - Inizializzati a zero
    std::vector<double> V(N, 0.0), V_new(N-1, 0.0);
    std::vector<double> I(N, 0.0), I_new(N-1, 0.0); // I[i] è tra i e i+1
    std::vector<double> J(N, 0.0); // input densità di corrente [A/m] dalla scarica degli rpc
    std::vector<double> time_factor(0.0);

    const bool savegif = false;
    if(savegif == true){
      // save profiles at t=0 for GIF animation
      save_profile(V, 0.0, N, dx, "rpc_signal");
      save_profile(I, 0.0, N-1, dx, "I");
      save_profile(J, 0.0, N, dx, "J"); // Passa N e dx
    }

    // Variabili per rilevamento bordi
    //const double threshold = std::abs(J_peak * Z0 * 0.1); // probably an error, (A/m)*ohm)=V/m cannot compare with V Soglia dinamica (es. 10% ampiezza stimata)
    const double threshold = aThreshold;//0.001; // V
    bool reached_left = false, reached_right = false;
    double time_left = -1.0, time_right = -1.0; // Inizializza a -1
    bool reached_left_two = false, reached_right_two = false;
    double time_left_two = -1.0, time_right_two = -1.0; // Inizializza a -1

    std::cout << "--- Threshold to detect signal at borders ---" << std::endl;
    std::cout << " threshold : " << threshold << " V " << std::endl;
    std::cout << "----------------------------------" << std::endl;

    // Time smearing
    double signaljitter = aSignaljitter;
    double TDCbinsize = aTDCbinsize;
    std::cout << "--- TDC resolution and time jitter ---" << std::endl;
    std::cout << " Signal time jitter : " << signaljitter << " s " 
              << " TDC bin size " << TDCbinsize << " s " << std::endl;
    std::cout << "----------------------------------" << std::endl;

    // Costanti moltiplicative per aggiornamento FDTD
    const double const_I = dt / (L * dx); // dt/(L*dx)
    const double const_V = dt / (C * dx); // dt/(C*dx)
    const double const_Source = dt / C;   // dt/C per termine sorgente

    #ifdef wLOSS
    const double denom_I = 1.0 + R * dt / (2.0 * L);
    const double cI1 = (1.0 - R * dt / (2.0 * L)) / denom_I;
    const double cI2 = (dt / (L * dx)) / denom_I;
    #endif

    #ifndef MURBC
    // <<< NUOVO BLOCCO: Costanti per BC Resistiva (G=0) >>>
    double cV1_bc = 0.0, cV2_bc = 0.0;
    const double R_load = R_L; // Usa R_L definita sopra
    if (std::abs(R_load) > 1e-18) { // Evita divisione per zero se R_L=0 (cortocircuito)
        const double c1_bc = C * dx / dt;
        const double c2_bc = 1.0 / (2.0 * R_load);
        const double den_bc = c1_bc + c2_bc;
        if (std::abs(den_bc) < 1e-18) {
             std::cerr << "Errore: den_bc (BC coeff) vicino a zero." << std::endl;
             std::abort();
        }
        cV1_bc = (c1_bc - c2_bc) / den_bc;
        cV2_bc = 1.0 / den_bc;
    } else {
        // Se R_L == 0 (cortocircuito), V al bordo è forzato a 0
        cV1_bc = 0.0;
        cV2_bc = 0.0; // V sarà 0
    }
    #endif

    // --- Parametri Emulazione Charge Sensitive Amplifier ---
    const double C_f_CSA = 0.2e-12;  // Capacità di feedback del CSA [F], ottenuto sapendo che il guadagno è 5mV/fC
    const double tau_CSA = 5.0e-9;  // Costante di tempo di decadimento CSA [s], ottenuto a BB5
    std::cout << "--- Parametri Charge Sensitive Amplifier Emulator ---" << std::endl;
    std::cout << " C_f_CSA = " << C_f_CSA * 1e12 << " pF, tau_CSA = " << tau_CSA * 1e9 << " ns" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    // Variabili per l'uscita del CSA emulato
    double V_CSA_left = 0.0;  // Uscita del CSA al bordo sinistro
    double V_CSA_right = 0.0; // Uscita del CSA al bordo destro
    // Vettori per salvare l'uscita del CSA nel tempo (opzionale, per plotting)
    std::vector<double> V_CSA_left_history;
    std::vector<double> V_CSA_right_history;
    std::vector<double> time_history_csa; // Per i tempi corrispondenti
    // Coefficienti per l'aggiornamento del CSA
    const double den_CSA = 1.0 + dt / (2.0 * tau_CSA);
    const double k_CSA1 = (1.0 - dt / (2.0 * tau_CSA)) / den_CSA;
    const double k_CSA2 = (dt / C_f_CSA) / den_CSA;
    
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
            int Jsteplimit = 0;
            for (int i = 1; i < N - 1; ++i) {
                double x = i * dx;
                double spatial_factor = std::exp(-std::pow((x - x_source), 2) / (2 * sigma_x * sigma_x));
                Jsteplimit = savegif ? steps : 400;
                if(step<=Jsteplimit && step % 20 ==0){ 
                  J[i] = J_peak * spatial_factor * temporal_factor;
                }
                double source_current_density = J_peak * spatial_factor * temporal_factor;

                // Applica a V: V_new[i] -= (dt / C) * J(x,t)
                V_new[i] -= const_Source * source_current_density;
            }
            if(step<=Jsteplimit && step % 20 ==0) save_profile(J, (step+0.5)*dt, N, dx, "J"); // Passa N e dx
        }

        #ifdef MURBC
        // Applica BC assorbenti di Mur per V[0] e V[N-1] (sovrascrive V_new ai bordi)
        // Usa V (valori a n) e V_new (valori a n+1 dei punti interni già calcolati)
        V_new[0] = V[1] + mur_const * (V_new[1] - V[0]);

        if (N > 1) {
            V_new[N - 1] = V[N - 2] + mur_const * (V_new[N - 2] - V[N - 1]);
        } else {
             V_new[0] = 0; // Caso N=1
        } // fine boundary condition
        #endif
 
        #ifndef MURBC
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

        // --- Emulazione Uscita CSA ---
        // Corrente che entra nel CSA al bordo sinistro (x=0)
        // I[0] è I_{1/2}, positiva se va a destra. Quindi la corrente IN un CSA a sx è -I[0]
        double I_strip_left = (N > 0) ? -I_new[0] : 0.0;

        // Corrente che entra nel CSA al bordo destro (x=L, i=N-1)
        // I[N-2] è I_{N-1/2}, positiva se va a destra. Questa è la corrente IN un CSA a dx.
        double I_strip_right = (N > 1) ? I_new[N-2] : 0.0;

        V_CSA_left = k_CSA1 * V_CSA_left - k_CSA2 * I_strip_left;
        V_CSA_right = k_CSA1 * V_CSA_right - k_CSA2 * I_strip_right;

        // Salva valori del CSA per plotting (opzionale)
        // t_output è il tempo alla fine dello step corrente ((step + 1) * dt)
        V_CSA_left_history.push_back(V_CSA_left);
        V_CSA_right_history.push_back(V_CSA_right);
        time_history_csa.push_back(t_output); // Assicurati che t_output sia definito qui
        // --- Fine Emulazione Uscita CSA ---

        // Controlla arrivo ai bordi (usa V_new appena calcolato)
        if (!reached_left && std::abs(V_CSA_left) >= threshold) { // Usa abs per soglia
            reached_left = true;
            time_left = t_output;
            // std::cout << "-> Segnale raggiunto bordo sinistro a t = " << t_output*1e9 << " ns" << std::endl;
        }
        if (!reached_right && std::abs(V_CSA_right) >= threshold) { // Usa abs per soglia
            reached_right = true;
            time_right = t_output;
            // std::cout << "-> Segnale raggiunto bordo destro a t = " << t_output*1e9 << " ns" << std::endl;
        } // fine controllo primo arrivo ai bordi
        
        // Controlla secondo arrivo ai bordi for time-over-threshold
        if (reached_left && !reached_left_two && std::abs(V_CSA_left) <= threshold) {
            reached_left_two = true;
            time_left_two = t_output;
        }
        if (reached_right && !reached_right_two && std::abs(V_CSA_right) <= threshold) {
            reached_right_two = true;
            time_right_two = t_output;
        }

        // Scrivi su file snapshot
        if (t_output < 1.5*tau && step % 20 == 0){ 
           save_profile(V_new, t_output, N, dx, "rpc_signal"); // Passa N e dx
           save_profile(I_new, t_output_forcurrent, N-1, dx, "I"); // Passa N e dx
        }
        if (savegif == true) snapshot_interval = 20;
        if ((t_output > 1.5*tau && step % snapshot_interval == 0) || step == steps -1 ) { // Salva anche l'ultimo
           // std::cout << "Salvataggio snapshot a t = " << t_output*1e9 << " ns (step " << step << ")" << std::endl;
           save_profile(V_new, t_output, N, dx, "rpc_signal"); // Passa N e dx
           save_profile(I_new, t_output_forcurrent, N-1, dx, "I"); // Passa N e dx
        }

        // Aggiorna stato per il prossimo ciclo
        V.swap(V_new);
        I.swap(I_new);

        // Condizione di uscita anticipata (opzionale)
        if (reached_left && reached_right && step > steps / 1 ) { // for the moment let's include all the steps step> steps / 1 (was / 2)
             double max_V_abs = 0.0;
             for(double val : V) { max_V_abs = std::max(max_V_abs, std::abs(val)); }
             if (max_V_abs < output_precision) {
                  std::cout << "\nSegnale uscito e tensione residua bassa a t = " << t_output*1e9 << " ns. Interrompo." << std::endl;
                  // Salva ultimo stato se non già salvato
                   if (step % snapshot_interval != 0) {
                       save_profile(V, t_output, N, dx, "rpc_signal");
                       save_profile(I, t_output_forcurrent, N-1, dx, "I");
                   }
                  break;
             }
        }

    } // Fine loop temporale

    save_timefactor(time_factor, dt);
    save_csa_output(time_history_csa, V_CSA_left_history, V_CSA_right_history);
    
    std::cout << "\n--- Risultati Rilevamento Bordi ---" << std::endl;
    if (time_left > 0)
       std::cout << "Tempo di arrivo a sinistra: " << time_left * 1e9 << " ns\n";
    else
       std::cout << "Segnale non rilevato a sinistra (soglia=" << threshold <<").\n";
    if (time_left_two > 0)
        std::cout << "Tempo di discesa a sinistra: " << time_left_two * 1e9 << "ns\n";
        std::cout << "Time over threshold a sinistra: " << (time_left_two - time_left) * 1e9 << "ns\n";

    if (time_right > 0)
       std::cout << "Tempo di arrivo a destra: " << time_right * 1e9 << " ns\n";
    else
       std::cout << "Segnale non rilevato a destra (soglia=" << threshold <<").\n";
    if (time_right_two > 0)
        std::cout << "Tempo di discesa a destra: " << time_right_two * 1e9 << "ns\n";
        std::cout << "Time over threshold a destra: " << (time_right_two - time_right) * 1e9 << "ns\n";

    std::cout << "Velocità teorica di propagazione: " << v << " m/s\n";

    rpcoutput thisOutput;
    thisOutput.add_metadata(length, N, R, tau, threshold, J_peak, signaljitter, TDCbinsize);
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<double> dist(0., signaljitter);
    // If signaljitter > 0 add it to time
    if (signaljitter > 0.){
        time_left = time_left + dist(gen);
        time_right = time_right + dist(gen);
        time_left_two = time_left_two + dist(gen);
        time_right_two = time_right_two + dist(gen);
    }
    std::uniform_real_distribution<double> uniform_dist(-TDCbinsize/2., TDCbinsize/2.);
    if(TDCbinsize > 0.){
        time_left = (std::round(time_left / TDCbinsize) * TDCbinsize) + uniform_dist(gen);
        time_right = (std::round(time_right / TDCbinsize) * TDCbinsize) + uniform_dist(gen);
        time_left_two = (std::round(time_left_two / TDCbinsize) * TDCbinsize) + uniform_dist(gen);
        time_right_two = (std::round(time_right_two / TDCbinsize) * TDCbinsize) + uniform_dist(gen);
    }
    thisOutput.add_output(time_left, time_right, time_left_two-time_left, time_right_two-time_right);
    thisOutput.set_filename(output_name);
    thisOutput.rpcoutput_tofile(append_output);

}

void parametrized_rpc(double aX, double aLength, double aChargeFraction, double signaljitter, double TDCbinsize, int NumberOfEvents){

    std::cout<<"Using parametrized rpc"<<std::endl;

    // this method uses parameterized TOA and TOT values,
    // add the effect of the signal jitter and TDC resolution
    // in order to study the effect on the phi resolution
    //
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<double> dist(0., signaljitter);
    std::uniform_real_distribution<double> uniform_dist(-TDCbinsize/2., TDCbinsize/2.);

    const double v = 0.23e9; // m/s
    const double dist_left = aX;
    const double dist_right = aLength-aX;
    double X_meas = 0.0;

    std::ostringstream filename;
    filename << "secondcoordinate" << std::to_string(aX) << "_" <<std::to_string(signaljitter)<<"_"<<std::to_string(TDCbinsize)<<".txt";
    std::ofstream out(filename.str());
     if (!out) {
        std::cerr << "Error while opening the file: " << filename.str() << std::endl;
        return;
    }
    out << std::fixed << std::setprecision(6);
    out <<" X_real (m) "<<" X_meas (m) "<< " jitter (ns) "<< " TDC size (ns) "<<std::endl;


    for(std::size_t i=0;i<NumberOfEvents; i++){
      
        double TOA_left = 5.9966*dist_left+0.3369*aChargeFraction-0.4742*std::pow(dist_left,2)-0.7377*dist_left*aChargeFraction - 0.2242*std::pow(aChargeFraction,2);
        double TOA_right = 5.9966*dist_right+0.3369*aChargeFraction-0.4742*std::pow(dist_right,2)-0.7377*dist_right*aChargeFraction - 0.2242*std::pow(aChargeFraction,2);
        if (signaljitter > 0.){
           TOA_left = TOA_left + dist(gen);
           TOA_right = TOA_right + dist(gen);
        }
        if(TDCbinsize > 0.){
           TOA_left = (std::round(TOA_left / TDCbinsize) * TDCbinsize) + uniform_dist(gen);
           TOA_right = (std::round(TOA_right / TDCbinsize) * TDCbinsize) + uniform_dist(gen);
        }
        X_meas = 1. + ((TOA_left/1e9-TOA_right/1e9)*v)/2. + aLength/v; //m
        std::cout<<"aX "<<aX<<" m "<<X_meas<<" m "<<std::endl;
        out << aX << "\t" << X_meas << "\t" << signaljitter << "\t" << TDCbinsize << "\n";
    }
}

void print_usage(const char* progName) {
    std::cout << "Uso: " << progName << " [--aLength valore (m)] [--aN valore] [--aR valore (ohm/m)] [--aTau valore (s)] [--aThreshold valore (V)] [--aJpeak valore (A/m)]\n"
              << "Tutti i parametri sono opzionali e hanno dei valori di default.\n";
}

int main(int argc, char* argv[]) {
    // Valori di default
    double aLength = 2.0; // m
    int aN = 4000;
    double aR = 0.02; // Ohm/m
    double aTau = 0.5e-9; // s
    double aThreshold = 0.01; // V
    double aJpeak = -7e-3; // A/m
    double signaljitter = 0.0;//8e-9; // s
    double TDCbinsize = 0.0;//8e-9; // s

    // Parsing degli argomenti
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--aLength" && i + 1 < argc) {
            aLength = atof(argv[++i]);
        } else if (arg == "--aN" && i + 1 < argc) {
            aN = atoi(argv[++i]);
        } else if (arg == "--aR" && i + 1 < argc) {
            aR = atof(argv[++i]);
        } else if (arg == "--aTau" && i + 1 < argc) {
            aTau = atof(argv[++i]);
        } else if (arg == "--aThreshold" && i + 1 < argc) {
            aThreshold = atof(argv[++i]);
        } else if (arg == "--aJpeak" && i + 1 < argc) {
            aJpeak = atof(argv[++i]);
        } else if (arg == "--aJitter" && i + 1 < argc) {
            signaljitter = atof(argv[++i]);
        } else if (arg == "--aTDCbin" && i + 1 < argc) {
            TDCbinsize = atof(argv[++i]);
        } 
        else {
            std::cerr << "Argomento non riconosciuto o mancante valore: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

   //parametrized_rpc(1.5, aLength, 1.0, 0.08, 0.08, 1000);

   process_rpc_signal(aLength, aN, aR, aTau, aThreshold, aJpeak, signaljitter, TDCbinsize);

    std::string outputname;    
    // Study behaviour as a function of threshold
    /*outputname = "threshold.txt";
    for(std::size_t i=0; i<10; i++){
        double newThreshold = 0.001 + i*0.0005; // V
        if(i==0) process_rpc_signal(aLength, aN, aR, aTau, newThreshold, aJpeak, signaljitter, TDCbinsize, outputname);
        else process_rpc_signal(aLength, aN, aR, aTau, newThreshold, aJpeak, signaljitter, TDCbinsize, outputname, true);
    }*/

    // Study behaviour as a function of tau
    /*outputname = "tau.txt";
    for(std::size_t i=0; i<30; i++){
        double newTau = 0.1e-9 + i*0.05e-9; // s
        if(i==0) process_rpc_signal(aLength, aN, aR, newTau, aThreshold, aJpeak, signaljitter, TDCbinsize, outputname);
        else process_rpc_signal(aLength, aN, aR, newTau, aThreshold, aJpeak, signaljitter, TDCbinsize, outputname, true);
    }*/

    // Study behaviour as a function of length
    /*outputname = "length.txt";
    for(std::size_t i=0; i<40; i++){
        double newLength = 1.0 + i*0.1; // m
        if(i==0) process_rpc_signal(newLength, aN, aR, aTau, aThreshold, aJpeak, signaljitter, TDCbinsize, outputname);
        else process_rpc_signal(newLength, aN, aR, aTau, aThreshold, aJpeak, signaljitter, TDCbinsize, outputname, true);
    }*/
    
    // Study behaviour as a function of Jpeak
    /*outputname = "jpeak.txt";
    for(std::size_t i=0; i<16; i++){
        double newJpeak = -0.001 - i*0.0004; // m
        if(i==0) process_rpc_signal(aLength, aN, aR, aTau, aThreshold, newJpeak, signaljitter, TDCbinsize, outputname);
        else process_rpc_signal(aLength, aN, aR, aTau, aThreshold, newJpeak, signaljitter, TDCbinsize, outputname, true);
    }*/

    // Study behaviour as a function ok Jpeak and Length
    /*outputname = "jpeak_length.txt";
    for(std::size_t i=0; i<15; i++){
        double newJpeak = -0.001 - i*0.0004; // m
        for(std::size_t j=0; j<40; j++){
            double newLength = 0.1 + j*0.05; // m
            if(i==0 && j==0) process_rpc_signal(newLength, aN, aR, aTau, aThreshold, newJpeak, signaljitter, TDCbinsize, outputname);
            else process_rpc_signal(newLength, aN, aR, aTau, aThreshold, newJpeak, signaljitter, TDCbinsize, outputname, true);
        }
    }*/

    return 0;
}
