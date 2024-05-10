#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <time.h>

using namespace std;

// Random number generator
random_device rd;
mt19937 rnd_gen;
uniform_real_distribution<double> uniDist(0.0,1.0);
normal_distribution<double> whiteNoise(0.0,1.0);

// Global variables
double ell, kappa, delta, x0, dt, dt_save, t, total_time;
int samples, seed;
unsigned long long steps_total, steps_save;

// Function declarations
void initialize(int argc, char **argv);
void save_position(double current_t, double pos, char phase);

// Main class of functions acting on the particle
class Particle {
public:
    double x;
    char phase;
    double t;
    double total_energy;

    Particle() {
        x = x0;
        phase = 'd';
        t = 0;
        total_energy = 0;
    }
    
    void pbc_wrap() {
        x = x - ell*floor(x/ell);
    }
    
    double get_force() {
        if (x < ell*delta) {
            return -kappa/(delta*ell);
        } else if (x > ell*delta) {
            return kappa/(ell*(1-delta));
        } else {
            return 0;
        }
    }

    double get_energy() {
        if (x <= ell*delta) {
            return kappa/(delta*ell) * x;
        } else if (x > ell*delta) {
            return kappa/(ell*(1-delta)) * (ell - x);
        } else {
            return 0;
        }
    }
    
    void move() {
        if (phase == 'd') {
            if (uniDist(rnd_gen) < dt) {
                phase = 'r';
                total_energy += get_energy();
                x += get_force() * dt + sqrt(2 * dt) * whiteNoise(rnd_gen);
            } else {
                x += sqrt(2 * dt) * whiteNoise(rnd_gen);
            }
            pbc_wrap();
        } else if (phase == 'r') {
            x += get_force() * dt + sqrt(2 * dt) * whiteNoise(rnd_gen);
            
            if (x <= 0 || x >= ell) {
                phase = 'd';
                pbc_wrap();
            }
        }
        else {
            cerr << "Error: particle phase is not 'd' or 'r'." << endl;
            exit(1);
        }

        t += dt;
    }
};

// Initialize parameters from input file
void initialize(int argc, char **argv) {
    std::fstream input_file;

    if (argc < 2) {
    	cerr << "Usage: " << argv[0] << " INPUT FILENAME" << endl;
    	exit(1);
	} else {
	    input_file.open(argv[1],ios::in);
	    if (input_file.fail()) 
	    {cerr << "Can't open input parameters file!" << endl; exit(1);}
	}

    input_file >> ell;
    input_file >> kappa;
    input_file >> delta;
    input_file >> x0;
    input_file >> dt;
    input_file >> dt_save;
    input_file >> total_time;
    input_file >> samples;
    input_file >> seed;
    input_file.close();

    steps_total = total_time / dt;
    steps_save = dt_save / dt;

    cout << "Initializing parameters..." << endl;
    cout << "ell = " << ell << ", kappa = " << kappa << ", delta = " << delta << endl;
    cout << "x0 = " << x0 << endl;
    cout << "dt = " << dt << ", total_time = " << total_time << ", dt_save = " << dt_save << endl;
    cout << "total steps = " << steps_total << ", save every " << steps_save << " steps" << endl;
    cout << "seed = " << seed << ", samples = " << samples << endl;

}

void save_position(double current_t, double pos, char phase) {
    std::fstream output_file;
    std::string p = string(1,phase);
    output_file.open(p + "_pos_" + to_string(current_t), ios::app);
    output_file << pos << endl;
    output_file.close();
}

void save_energy(double current_t, double total_energy) {
    std::fstream output_file;
    output_file.precision(19);
    output_file.open("energy_" + to_string(current_t), ios::app);
    output_file << total_energy << endl;
    output_file.close();
}

int main(int argc, char **argv) {

    initialize(argc, argv);

    // Open output files
    for (double i = dt_save; i <= total_time; i += dt_save) {
        std::fstream output_file;
        output_file.precision(19);
        output_file.open("energy_" + to_string(i), ios::out);
        if (output_file.fail())
        {cerr << "Can't open output file!" << endl; exit(1);}
        output_file.close();
    }

    rnd_gen.seed (seed);

    // Main simulation loop
    double current_t = 0;
    cout << "Starting simulation ..." << endl;
    cout << "Progress: " << flush;
    for (int n = 0; n < samples; n++) {
        Particle particle;
        for (unsigned long long i = 0; i <= steps_total; i++) {
            if (i % steps_save == 0 && i != 0) {
                current_t = i * dt;
                current_t = round(current_t);
                save_energy(current_t, particle.total_energy);
            }
            particle.move();
        }
        if ((n+1) % int(floor(samples / 10)) == 0) {
            cout << "|" << flush;
        }
    }

    

    return 0;
}

