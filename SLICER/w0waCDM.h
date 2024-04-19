/**
 * @file w0waCDM.h
 * @brief Class to model the universe using the w0waCDM cosmological model.
 *
 * This class provides methods to compute cosmological distances such as
 * comoving distance and transverse comoving distance under the w0waCDM model,
 * which includes parameters for the equation of state of dark energy that
 * vary with redshift.
 */

#ifndef W0WACDM_H
#define W0WACDM_H

#include <cmath>
#include <map>
#include "utilities.h"

/**
 * @class w0waCDM
 * @brief A class to compute distances in a universe modeled by the w0waCDM cosmology.
 *
 * This class encapsulates the parameters and functions required to compute
 * the Hubble parameter and various cosmological distances such as comoving
 * and transverse comoving distances considering a dynamic dark energy component.
 */
class w0waCDM {
private:
    static constexpr double CSPEEDOFLIGHT = speedcunit * 100; ///< Speed of light in H0/h * 1 Mpc
    double H0;  ///< Hubble constant at z = 0 in units of km/s/Mpc
    double OmegaM;  ///< Matter density parameter
    double OmegaLambda;  ///< Dark energy density parameter
    double w0;  ///< Equation of state parameter at z = 0
    double wa;  ///< Change rate of the equation of state parameter
    mutable std::map<double, double> cache;  ///< Cache for storing computed comoving distances

    /**
     * Calculates the Hubble parameter H(z) as a function of redshift z.
     * @param z The redshift at which to compute the Hubble parameter.
     * @return The Hubble parameter in km/s/Mpc.
     */
    double Hz(double z) const;

public:
    /**
     * Constructs a w0waCDM cosmology model with specified parameters.
     * @param H0 Hubble constant at z = 0 in km/s/Mpc.
     * @param OmegaM Matter density parameter.
     * @param OmegaLambda Dark energy density parameter.
     * @param w0 Equation of state parameter at z = 0.
     * @param wa Change rate of the equation of state parameter.
     */
    w0waCDM(double H0, double OmegaM, double OmegaLambda, double w0, double wa);

    /**
     * Computes the comoving distance from z = 0 to a given redshift z.
     * This method uses caching to speed up repeated calculations.
     * @param z The redshift to which the comoving distance is calculated.
     * @return The comoving distance in Mpc.
     */
    double comovingDistance(double z) const;

    /**
     * Computes the transverse comoving distance at a given redshift z.
     * This distance accounts for the curvature of space.
     * @param z The redshift at which to compute the transverse comoving distance.
     * @return The transverse comoving distance in Mpc.
     */
    double transverseComovingDistance(double z) const;
};

#endif // W0WACDM_H