#ifndef INCLUDE_VISCOUS_FLUX_H_
#define INCLUDE_VISCOUS_FLUX_H_

#include "GasModel.h"

/* no_viscous_flux
 * functor that computes no viscous flux for inviscid calculation
 * Uses the enum isViscous=false to remove unneeded code at compile time.
 */
struct no_viscous_flux {
  enum { isViscous = false };
  typedef typename ViewTypes::c_rnd_vector_field_type vector_field_type;

  no_viscous_flux() {}

  void compute_flux(const double grad_primitive[5][3],
                    const double* const primitive, const double* const& a_vec,
                    double* vflux) const {
    return;
  }
};

/* newtonian_viscous_flux
 * functor that computer Newtonian viscous flux using gradients
 * and primitive values.
 */
struct newtonian_viscous_flux {
  enum { isViscous = true };

  newtonian_viscous_flux() {}

  void compute_flux(const double grad_primitive[5][3],
                    const double* const& primitive, const double* const& a_vec,
                    double* const& vflux) const {
    double viscosity = ComputeViscosity(primitive[4]);
    double thermal_conductivity = ComputeThermalConductivity(viscosity);
    double divergence_velocity = 0;

    for (int icomp = 0; icomp < 5; ++icomp) {
      vflux[icomp] = 0.0;
    }

    for (int i = 0; i < 3; ++i) {
      divergence_velocity += grad_primitive[i + 1][i];
    }

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const double delta_ij = (i == j) ? 1 : 0;
        const double S_ij =
            0.5 * (grad_primitive[i + 1][j] + grad_primitive[j + 1][i]);
        const double t_ij = S_ij - divergence_velocity * delta_ij / 3.;
        vflux[1 + i] += (2 * viscosity * t_ij) * a_vec[j];
        vflux[4] += (2 * viscosity * t_ij) * primitive[i + 1] * a_vec[j];
      }

      vflux[4] += thermal_conductivity * grad_primitive[4][i] * a_vec[i];
    }

    return;
  }
};

#endif
