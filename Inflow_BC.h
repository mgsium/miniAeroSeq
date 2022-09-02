#ifndef INCLUDE_INFLOW_BC_H_
#define INCLUDE_INFLOW_BC_H_

#include <cstdio>
#include "Faces.h"

/*compute_inflow_BC_flux
 * functor to compute contribution of inflow boundary condition
 * state is fixed and passed as flow_state
 * The state along with the internal cell state are used to compute
 * the external face flux contribution
 */
template<class FluxType>
struct compute_inflowBC_flux {
  typedef typename ViewTypes::c_rnd_solution_field_type solution_field_type;
  typedef typename ViewTypes::c_rnd_face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  typedef typename ViewTypes::c_vector_field_type vector_field_type;

  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  cell_storage_field_type cell_flux_;
  vector_field_type face_normal_, face_tangent_, face_binormal_;
  FluxType flux_evaluator_;
  double flow_state_[5];

  compute_inflowBC_flux(Faces faces, solution_field_type cell_values,
      Cells cells, double * flow_state, FluxType flux) :
      face_cell_conn_(faces.face_cell_conn_), cell_flux_index_(
          faces.cell_flux_index_), cell_values_(cell_values), cell_flux_(
          cells.cell_flux_), face_normal_(faces.face_normal_), face_tangent_(
          faces.face_tangent_), face_binormal_(faces.face_binormal_), flux_evaluator_(
          flux) {
    for (int i = 0; i < 5; ++i)
      flow_state_[i] = flow_state[i];
  }

  void operator()(int i) const {
    int index = face_cell_conn_[i][0];

    double flux[5];
    double conservatives_r[5];
    double conservatives_l[5];
    double primitives_r[5];
    double primitives_l[5];

    for (int icomp = 0; icomp < 5; ++icomp) {
      conservatives_l[icomp] = cell_values_[index][icomp];
      conservatives_r[icomp] = flow_state_[icomp];
    }

    ComputePrimitives(conservatives_l, primitives_l);
    ComputePrimitives(conservatives_r, primitives_r);


    flux_evaluator_.compute_flux(primitives_l, primitives_r, flux, &face_normal_[i][0],
        &face_tangent_[i][0], &face_binormal_[i][0]);

#ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; ++icomp)
    {
      cell_flux_[index][cell_flux_index_[i][0]][icomp] = -flux[icomp];
    }
#endif
  }
};

#endif