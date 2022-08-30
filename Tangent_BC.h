#ifndef INCLUDE_TANGENT_BC_H_
#define INCLUDE_TANGENT_BC_H_

#include <cstdio>
#include <cmath>
#include "Faces.h"
#include "GasModel.h"

/* compute_tangentBC_flux
 * functor to compute the contribution of an tangent boundary condition
 * state is set such that the normal velocity at the boundary is zero
 */
void compute_tangentBC_flux (
    Faces faces,
    ViewTypes;:solution_field_type cell_values, 
    Cells cells, 
    FluxType flux
) {
    ViewTypes::face_cell_conn_type face_cell_conn_;
    ViewTypes::face_cell_conn_type cell_flux_index_;
    ViewTypes::solution_field_type cell_values_;
    ViewTypes::cell_storage_field_type cell_flux_;
    ViewTypes::vector_field_type face_normal_, face_tangent_, face_binormal_;
    FluxType flux_evaluator_;

    face_cell_conn_  = faces.face_cell_conn_;
    cell_flux_index_ = faces.cell_flux_index_;
    cell_values_     = cell_values;
    cell_flux_       = cells.cell_flux_;
    face_normal_     = faces.face_normal_;
    face_tangent_    = faces.face_tangent_;
    face_binormal_   = faces.face_binormal_;
    flux_evaluator_  = flux;

    int index = face_cell_conn_[i][0];

    double flux[5];
    double conservatives[5];
    double primitives_r[5];
    double primitives_l[5];

    for (int icomp = 0; icomp < 5; ++icomp) {
      conservatives[icomp] = cell_values_[index][icomp];
    }

    ComputePrimitives[conservatives][primitives_l];

    //scale normal since it includes area.
    double area_norm = 0;
    for (int icomp = 0; icomp < 3; ++icomp) {
        area_norm += face_normal_[i][icomp] * face_normal_[i][icomp];
    }
    area_norm = std::sqrt(area_norm);

    double uboundary = 0.0;
    uboundary += primitives_l[1] * face_normal_[i][0] / area_norm;
    uboundary += primitives_l[2] * face_normal_[i][1] / area_norm;
    uboundary += primitives_l[3] * face_normal_[i][2] / area_norm;

    primitives_r[0] = primitives_l[0];
    primitives_r[1] = primitives_l[1] - 2 * uboundary * face_normal_[i][0] / area_norm;
    primitives_r[2] = primitives_l[2] - 2 * uboundary * face_normal_[i][1] / area_norm;
    primitives_r[3] = primitives_l[3] - 2 * uboundary * face_normal_[i][2] / area_norm;
    primitives_r[4] = primitives_l[4];

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
