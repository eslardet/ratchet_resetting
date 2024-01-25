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
double L, a, h, D, r, x0, dt, dt_save, t, total_time;
int samples, steps_total, steps_save, seed;

// Function declarations
void initialize(int argc, char **argv);
void save_position(int current_t, float pos, char phase);

// Main class of functions acting on the particle
class Particle {
public:
    double x;
    char phase;
    double t;

    Particle() {
        x = x0;
        phase = 'd';
        t = 0;
    }
    
    void pbc_wrap() {
        x = x - L*floor(x/L);
    }
    
    double get_force() {
        if (x < a) {
            return -h / a;
        } else if (x > a) {
            return h / (L - a);
        } else {
            return 0;
        }
    }
    
    void move() {
        if (phase == 'd') {
            if (uniDist(rnd_gen) < r * dt) {
                phase = 'r';
                x += get_force() * dt + sqrt(2 * D * dt) * whiteNoise(rnd_gen);
            } else {
                x += sqrt(2 * D * dt) * whiteNoise(rnd_gen);
            }
            pbc_wrap();
        } else if (phase == 'r') {
            x += get_force() * dt + sqrt(2 * D * dt) * whiteNoise(rnd_gen);
            
            if (x <= 0 || x >= L) {
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

    input_file >> L;
    input_file >> a;
    input_file >> h;
    input_file >> D;
    input_file >> r;
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
    cout << "L = " << L << ", a = " << a << ", h = " << h << endl;
    cout << "D = " << D << ", r = " << r << ", x0 = " << x0 << endl;
    cout << "dt = " << dt << ", total_time = " << total_time << ", dt_save = " << dt_save << endl;
    cout << "total steps = " << steps_total << ", save every " << steps_save << " steps" << endl;
    cout << "seed = " << seed << ", samples = " << samples << endl;

}

void save_position(int current_t, float pos, char phase) {
    std::fstream output_file;
    std::string p = string(1,phase);
    output_file.open(p + "_pos_" + to_string(current_t), ios::app);
    output_file << pos << endl;
    output_file.close();
}

int main(int argc, char **argv) {

    initialize(argc, argv);

    // Open output files
    char phases[2] = {'d', 'r'};
    for (int i = 0; i <= total_time; i += dt_save) {
        for (int j = 0; j < 2; j++) {
            std::string p = string(1, phases[j]);
            std::fstream output_file;
            output_file.open(p + "_pos_" + to_string(i), ios::out);
            if (output_file.fail())
            {cerr << "Can't open output file!" << endl; exit(1);}
            output_file.close();
        }
    }

    rnd_gen.seed (seed);

    // Main simulation loop
    int current_t = 0;
    for (int n = 0; n < samples; n++) {
        Particle particle;
        for (int i = 0; i <= steps_total; i++) {
            if (i % steps_save == 0) {
                current_t = i * dt;
                save_position(current_t, particle.x, particle.phase);
            }
            particle.move();
        }
    }

    

    return 0;
}

