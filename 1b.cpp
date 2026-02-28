#include <iostream>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <iomanip>

using namespace std;

struct KeyHash {
    size_t operator()(const uint64_t& k) const noexcept { return std::hash<uint64_t>{}(k); }
};

// Pack (x1,x2,x3) into 64-bit key (supports up to 2^21-1 per component)
static inline uint64_t pack(int x1, int x2, int x3) {
    return ( (uint64_t)(x1 & 0x1FFFFF) << 42 ) |
           ( (uint64_t)(x2 & 0x1FFFFF) << 21 ) |
           ( (uint64_t)(x3 & 0x1FFFFF) );
}

static inline void unpack(uint64_t k, int &x1, int &x2, int &x3) {
    x1 = (int)((k >> 42) & 0x1FFFFF);
    x2 = (int)((k >> 21) & 0x1FFFFF);
    x3 = (int)( k        & 0x1FFFFF);
}

int main() {
    const int STEPS = 7;

    unordered_map<uint64_t, double, KeyHash> dist, nextDist;
    dist.reserve(512);
    nextDist.reserve(512);

    // Start: [9,8,7]
    dist[pack(9,8,7)] = 1.0;

    for (int step = 0; step < STEPS; step++) {
        nextDist.clear();

        for (const auto &kv : dist) {
            int x1,x2,x3;
            unpack(kv.first, x1,x2,x3);
            double p = kv.second;

            // propensities from the homework discrete model:
            // a1 = (1/2) x1(x1-1) x2
            // a2 = x1 x3(x3-1)
            // a3 = 3 x2 x3
            double a1 = 0.5 * (double)x1 * (x1 - 1) * (double)x2;
            double a2 = (double)x1 * (double)x3 * (x3 - 1);
            double a3 = 3.0 * (double)x2 * (double)x3;
            double sum = a1 + a2 + a3;

            if (sum == 0.0) {
                nextDist[kv.first] += p;
                continue;
            }

            // R1: 2X1 + X2 -> 4X3
            if (a1 > 0.0) {
                int nx1 = x1 - 2, nx2 = x2 - 1, nx3 = x3 + 4;
                if (nx1 >= 0 && nx2 >= 0 && nx3 >= 0)
                    nextDist[pack(nx1,nx2,nx3)] += p * (a1 / sum);
            }

            // R2: X1 + 2X3 -> 3X2
            if (a2 > 0.0) {
                int nx1 = x1 - 1, nx2 = x2 + 3, nx3 = x3 - 2;
                if (nx1 >= 0 && nx2 >= 0 && nx3 >= 0)
                    nextDist[pack(nx1,nx2,nx3)] += p * (a2 / sum);
            }

            // R3: X2 + X3 -> 2X1
            if (a3 > 0.0) {
                int nx1 = x1 + 2, nx2 = x2 - 1, nx3 = x3 - 1;
                if (nx1 >= 0 && nx2 >= 0 && nx3 >= 0)
                    nextDist[pack(nx1,nx2,nx3)] += p * (a3 / sum);
            }
        }

        dist.swap(nextDist);
    }

    // Compute means and variances from the exact distribution after 7 steps
    double EX1=0, EX2=0, EX3=0;
    double EX1_2=0, EX2_2=0, EX3_2=0;
    double totalProb=0;

    for (const auto &kv : dist) {
        int x1,x2,x3;
        unpack(kv.first, x1,x2,x3);
        double p = kv.second;

        totalProb += p;
        EX1 += x1 * p; EX2 += x2 * p; EX3 += x3 * p;
        EX1_2 += (double)x1 * x1 * p;
        EX2_2 += (double)x2 * x2 * p;
        EX3_2 += (double)x3 * x3 * p;
    }

    // Normalize (should already be ~1, but be robust)
    EX1 /= totalProb; EX2 /= totalProb; EX3 /= totalProb;
    EX1_2 /= totalProb; EX2_2 /= totalProb; EX3_2 /= totalProb;

    double Var1 = EX1_2 - EX1*EX1;
    double Var2 = EX2_2 - EX2*EX2;
    double Var3 = EX3_2 - EX3*EX3;

    cout << fixed << setprecision(6);
    cout << "After " << STEPS << " steps from [9,8,7] (exact):\n";
    cout << "E[X1] = " << EX1 << "   Var(X1) = " << Var1 << "\n";
    cout << "E[X2] = " << EX2 << "   Var(X2) = " << Var2 << "\n";
    cout << "E[X3] = " << EX3 << "   Var(X3) = " << Var3 << "\n";
    cout << "Total probability mass = " << totalProb << "\n";
    cout << "Reachable states counted = " << dist.size() << "\n";

    return 0;
}