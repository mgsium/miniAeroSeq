#ifndef INCLUDE_GASMODEL_H_
#define INCLUDE_GASMODEL_H_
#include <cmath>

/*Functions to compute gas properties from the conserved variables
 * or primitive variables using an ideal gas model.
 */
double ComputePressure(const double* V)
{
    const double Rgas = 287.05;
    const double rho = V[0];
    const double T = V[4];

    return rho*Rgas*T;
}

// ------------------------------------------------------------------------------------------------

double ComputeSoundSpeed(const double* V)
{
    const double gamma = 1.4;
    const double rho = V[0];

    const double pressure = ComputePressure(V);

    return std::sqrt(gamma * pressure / rho);
}

// ------------------------------------------------------------------------------------------------

double ComputeEnthalpy(const double* V)
{
    const double Cp = 1004.0;
    const double T = V[4];
    return Cp*T;
}

void ComputePrimitives(const double* U, double* V)
{
    double gamma = 1.4;
    double Rgas = 287.05;
    double r, u, v, w, T, ri, k, e;

    r  = U[0];
    ri = 1.0 / r;
    u  = U[1] * ri;
    v  = U[2] * ri;
    w  = U[3] * ri;
    k  = 0.5 * (u * u + v * v + w * w);
    e  = U[4] * ri - k;
    T  = e * (gamma - 1.0) / Rgas;

    V[0] = r;
    V[1] = u;
    V[2] = v;
    V[3] = w;
    V[4] = T;
}

double ComputeViscosity(const double temperature)
{
    const double sutherland_0 = 1.458e-6;
    const double sutherland_1 = 110.4;
    return sutherland_0 * temperature * std::sqrt(temperature) / (temperature + sutherland_1);
}

double ComputeThermalConductivity(const double viscosity) 
{
    const double Pr = 0.71;
    const double Cp = 1006.0;
    return viscosity*Cp/Pr; 
}

#endif
