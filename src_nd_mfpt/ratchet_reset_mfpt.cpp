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
double ell, kappa, delta, x0, dt;
int samples, seed;

// Function declarations
void initialize(int argc, char **argv);
void save_time(double current_t, double pos, char phase);

// Main class of functions acting on the particle
class Particle {
public:
    double x;
    char phase;
    double t;

    Particle() {
        x = x0;
        phase = 'r';
        t = 0;
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
    
    void move() {
        x += get_force() * dt + sqrt(2 * dt) * whiteNoise(rnd_gen);
        
        if (x <= 0 || x >= ell) {
            phase = 'd';
        }
        else {
            t += dt;
        }
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
    input_file >> samples;
    input_file >> seed;
    input_file.close();

    cout << "Initializing parameters..." << endl;
    cout << "ell = " << ell << ", kappa = " << kappa << ", delta = " << delta << endl;
    cout << "x0 = " << x0 << endl;
    cout << "dt = " << dt << endl;
    cout << "seed = " << seed << ", samples = " << samples << endl;

}

void save_time(double t_mfpt) {
    std::fstream output_file;
    output_file.open("mfpt_" + to_string(x0), ios::app);
    output_file << t_mfpt << endl;
    output_file.close();
}

int main(int argc, char **argv) {

    initialize(argc, argv);

    // Open output files
    std::fstream output_file;
    output_file.open("mfpt_" + to_string(x0), ios::out);
            if (output_file.fail())
            {cerr << "Can't open output file!" << endl; exit(1);}
            output_file.close();

    rnd_gen.seed (seed);

    // Main simulation loop
    cout << "Starting simulation ..." << endl;
    cout << "Progress: " << flush;
    for (int n = 0; n < samples; n++) {
        Particle particle;
        while (particle.phase == 'r') {
            particle.move();
        }
        save_time(particle.t);

        if ((n+1) % int(floor(samples / 10)) == 0) {
            cout << "|" << flush;
        }
    }

    

    return 0;
}

