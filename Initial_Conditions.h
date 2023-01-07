#ifndef INCLUDE_INITIAL_CONDITIONS_H_
#define INCLUDE_INITIAL_CONDITIONS_H_
#include "Cells.h"

/* initialize_sod3d
 * functor to set initialize conditions
 * sets up the SOD problem (shock-tube problem)
 * which is really 1D but can be run 2D or 3D.
 */
void initialize_sod3d(Cells cells_, ViewTypes::solution_field_type soln_,
                      ViewTypes::solution_field_type solnp1_, double midx_,
                      int i) {
  const double Rgas = 287.05;
  const double gamma = 1.4;
  const double Cv = Rgas / (gamma - 1.0);

  double P1 = 68947.57;
  double T1 = 288.889;
  double P2 = 6894.757;
  double T2 = 231.11;

  double density1 = P1 / (Rgas * T1);
  double rhoE1 = density1 * (Cv * T1);

  double density2 = P2 / (Rgas * T2);
  double rhoE2 = density2 * (Cv * T2);
  double x = cells_.coordinates_[0][i];

  if (x < midx_) {
    soln_[0][i] = density1;
    soln_[1][i] = 0.0;
    soln_[2][i] = 0.0;
    soln_[3][i] = 0.0;
    soln_[4][i] = rhoE1;
    for (int icomp = 0; icomp < 5; icomp++) {
      solnp1_[icomp][i] = soln_[icomp][i];
    }
  } else {
    soln_[0][i] = density2;
    soln_[1][i] = 0.0;
    soln_[2][i] = 0.0;
    soln_[3][i] = 0.0;
    soln_[4][i] = rhoE2;
    for (int icomp = 0; icomp < 5; icomp++) {
      solnp1_[icomp][i] = soln_[icomp][i];
    }
  }
}

/* initialize_constant
 * functor to set initialize conditions
 * set constant condition for entire flow field
 * This constant state is passed in as flow_state
 */
struct initialize_constant {
  typedef typename ViewTypes::solution_field_type solution_field_type;

  struct Cells cells_;
  // solution_field_type soln_, solnp1_;
  double flow_state_[5];

  initialize_constant(struct Cells cells, double* flow_state) : cells_(cells) {
    for (int i = 0; i < 5; ++i) {
      flow_state_[i] = flow_state[i];
    }
  }

  void operator()(int i, solution_field_type soln_,
                  solution_field_type solnp1_) const {
    for (int j = 0; j < 5; j++) {
      soln_[j][i] = flow_state_[j];
    }

    for (int icomp = 0; icomp < 5; icomp++) {
      solnp1_[icomp][i] = soln_[icomp][i];
    }
  }
};

#endif