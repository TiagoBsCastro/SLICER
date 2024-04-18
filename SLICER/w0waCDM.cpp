#include <cmath>
#include <map>
#include <stdexcept>  // Include for std::invalid_argument
#include "w0waCDM.h"
#include "utilities.h"


// Constructor with parameter validation
w0waCDM::w0waCDM(double H0, double OmegaM, double OmegaLambda, double w0, double wa)
    : H0(H0), OmegaM(OmegaM), OmegaLambda(OmegaLambda), w0(w0), wa(wa) {
    if (H0 <= 0 || OmegaM < 0 || OmegaLambda < 0) {
        throw std::invalid_argument("Invalid cosmological parameters: H0 must be positive, and density parameters cannot be negative.");
    }
}

// Hz function implementation
double w0waCDM::Hz(double z) const {
    double rhoLambda = OmegaLambda * pow(1 + z, 3 * (1 + w0 + wa)) * exp(-3 * wa * z / (1 + z));
    double rhoM = OmegaM * pow(1 + z, 3);
    double rhoTot = rhoLambda + rhoM + (1 - OmegaM - OmegaLambda) * pow(1 + z, 2);
    return H0 * sqrt(rhoTot);
}

// Compute comoving distance from z = 0 to z = z
double w0waCDM::comovingDistance(double z) const {
    // Check if the value is already in the cache
    if (cache.find(z) != cache.end()) {
        return cache[z];
    }

    double distance = 0;
    double dz = 0.001;  // Integration step size
    double lastZ = 0;

    // Find the largest z in the cache less than the target z, if any
    auto it = cache.lower_bound(z);
    if (it != cache.begin()) {
        --it;
        distance = it->second;
        lastZ = it->first;
    }

    for (double zi = lastZ; zi < z; zi += dz) {
        distance += 0.5 * dz * (1.0 / Hz(zi) + 1.0 / Hz(zi + dz));
    }

    // Store the computed value in the cache
    cache[z] = distance * CSPEEDOFLIGHT;  // Convert to Mpc
    return cache[z];
}

    // Method to compute the transverse comoving distance D_M(z)
double w0waCDM::transverseComovingDistance(double z) const {
    double D_C = comovingDistance(z); // Radial comoving distance
    if (fabs(1 - OmegaM - OmegaLambda) < 1e-5) { // Flat universe
        return D_C;
    } else {
        double OmegaK = 1.0 - OmegaM - OmegaLambda;
        double sqrtOmegaK = sqrt(fabs(OmegaK));
        double D_M;

        if (OmegaK < 0) {
            D_M = CSPEEDOFLIGHT / H0 / sqrtOmegaK * sinh(sqrtOmegaK * H0 / CSPEEDOFLIGHT * D_C);
        } else {
            D_M = CSPEEDOFLIGHT / H0 / sqrtOmegaK * sin(sqrtOmegaK * H0 / CSPEEDOFLIGHT * D_C);
        }

        return D_M;
    }
}
