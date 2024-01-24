#include <iostream>
#include <random>
#include <cmath>
#include <time.h>
using namespace std;


class Particle {
public:
    double x;
    char phase;
    double t;
    double dt;
    double D;
    double r;
    double L;
    double a;
    double h;
    
    Particle(double x0, double dt, double D, double L, double a, double h, double r) {
        x = x0;
        phase = 'd';
        t = 0;
        dt = dt;
        D = D;
        r = r;
        L = L;
        a = a;
        h = h;
    }
    
    void pbc_wrap() {
        x = fmod(x, L);
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
            double prob_r = ((double) rand() / (RAND_MAX));
            if (prob_r < r * dt) {
                phase = 'r';
                x += get_force() * dt + sqrt(2 * D * dt) * ((double) rand() / (RAND_MAX));
            } else {
                x += sqrt(2 * D * dt) * ((double) rand() / (RAND_MAX));
            }
            pbc_wrap();
        } else if (phase == 'r') {
            x += get_force() * dt + sqrt(2 * D * dt) * ((double) rand() / (RAND_MAX));
            
            if (x <= 0 || x >= L) {
                phase = 'd';
                pbc_wrap();
            }
        }
        t += dt;
    }
};

int main() {
    clock_t tStart = clock();

    double L = 1;
    double a = 0.2;
    double h = 1;
    double D = 0.1;
    double dt = 0.001;
    double r = 0.2;
    double x0 = 0;

    int samples = pow(10,5);
    int total_time = 10;
    int steps = total_time / dt;

    for (int j = 0; j < samples; j++) {
        Particle particle(x0, dt, D, L, a, h, r);
        for (int i = 0; i < steps; i++) {
            particle.move();
        }

    }


    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    return 0;
}

