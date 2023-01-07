#ifndef INCLUDE_TANGENT_BC_H_
#define INCLUDE_TANGENT_BC_H_

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "Faces.h"
#include "GasModel.h"

/* compute_tangentBC_flux
 * functor to compute the contribution of an tangent boundary condition
 * state is set such that the normal velocity at the boundary is zero
 */
template <class FluxType>
struct compute_tangentBC_flux {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::c_rnd_face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  cell_storage_field_type cell_flux_;
  vector_field_type face_normal_, face_tangent_, face_binormal_;
  FluxType flux_evaluator_;

  compute_tangentBC_flux(Faces faces, solution_field_type cell_values,
                         Cells cells, FluxType flux)
      : flux_evaluator_(flux) {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2,
              face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2,
              cell_flux_index_);
    std::copy(cell_values, cell_values + 5, cell_values_);
    std::copy(cells.cell_flux_, cells.cell_flux_ + 5, cell_flux_);
    // std::copy(faces.face_normal_, faces.face_normal_ + 3, face_normal_);
    // std::copy(faces.face_tangent_, faces.face_tangent_ + 3, face_tangent_);
    // std::copy(faces.face_binormal_, faces.face_binormal_ + 3,
    // face_binormal_);
  }

  void operator()(int i, Faces faces) const {
    int index = face_cell_conn_[0][i];

    double flux[5];
    double conservatives[5];
    double primitives_r[5];
    double primitives_l[5];

    // printf("here %d\n", index);

    for (int icomp = 0; icomp < 5; icomp++) {
      conservatives[icomp] = cell_values_[icomp][index];
    }

    ComputePrimitives(conservatives, primitives_l);

    // scale normal since it includes area.
    double area_norm = 0;
    for (int icomp = 0; icomp < 3; icomp++) {
      area_norm += faces.face_normal_[icomp][i] * faces.face_normal_[icomp][i];
    }
    area_norm = std::sqrt(area_norm);

    double uboundary = 0.0;
    uboundary += primitives_l[1] * faces.face_normal_[0][i] / area_norm;
    uboundary += primitives_l[2] * faces.face_normal_[1][i] / area_norm;
    uboundary += primitives_l[3] * faces.face_normal_[2][i] / area_norm;

    primitives_r[0] = primitives_l[0];
    primitives_r[1] =
        primitives_l[1] - 2 * uboundary * faces.face_normal_[0][i] / area_norm;
    primitives_r[2] =
        primitives_l[2] - 2 * uboundary * faces.face_normal_[1][i] / area_norm;
    primitives_r[3] =
        primitives_l[3] - 2 * uboundary * faces.face_normal_[2][i] / area_norm;
    primitives_r[4] = primitives_l[4];

    flux_evaluator_.compute_flux(
        primitives_l, primitives_r, flux, &faces.face_normal_[0][i],
        &faces.face_tangent_[0][i], &faces.face_binormal_[0][i]);

#ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; icomp++) {
      cell_flux_[index][cell_flux_index_[0][i]][icomp] = -flux[icomp];
    }
#endif
  }
};

#endif
