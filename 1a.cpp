#include <iostream>
#include <random>
#include <vector>
#include <cmath>

using namespace std;

struct State {
    long x1;
    long x2;
    long x3;
};

int main() {
    const int TRIALS = 10000;   // increase for better accuracy
    const int MAX_STEPS = 100000;

    mt19937_64 rng(random_device{}());
    uniform_real_distribution<double> uni(0.0, 1.0);

    int countC1 = 0;
    int countC2 = 0;
    int countC3 = 0;

    for (int t = 0; t < TRIALS; t++) {

        State s = {110, 26, 55};

        for (int step = 0; step < MAX_STEPS; step++) {

            // Check stopping conditions
            if (s.x1 >= 150) { countC1++; break; }
            if (s.x2 < 10)   { countC2++; break; }
            if (s.x3 > 100)  { countC3++; break; }

            // Compute propensities (discrete model)
            double a1 = 0.5 * s.x1 * (s.x1 - 1) * s.x2;
            double a2 = s.x1 * s.x3 * (s.x3 - 1);
            double a3 = 3.0 * s.x2 * s.x3;

            double sum = a1 + a2 + a3;
            if (sum == 0) break;

            double r = uni(rng) * sum;

            // Select reaction
            if (r < a1) {
                // R1: 2X1 + X2 → 4X3
                s.x1 -= 2;
                s.x2 -= 1;
                s.x3 += 4;
            }
            else if (r < a1 + a2) {
                // R2: X1 + 2X3 → 3X2
                s.x1 -= 1;
                s.x3 -= 2;
                s.x2 += 3;
            }
            else {
                // R3: X2 + X3 → 2X1
                s.x2 -= 1;
                s.x3 -= 1;
                s.x1 += 2;
            }

            // prevent negative states
            if (s.x1 < 0 || s.x2 < 0 || s.x3 < 0)
                break;
        }
    }

    cout << "Estimated probabilities:\n";
    cout << "Pr(C1) = " << (double)countC1 / TRIALS << endl;
    cout << "Pr(C2) = " << (double)countC2 / TRIALS << endl;
    cout << "Pr(C3) = " << (double)countC3 / TRIALS << endl;

    return 0;
}