#ifndef INCLUDE_NOSLIP_BC_H_
#define INCLUDE_NOSLIP_BC_H_

#include <cstdio>
#include <cmath>
#include "Faces.h"
#include "GasModel.h"


/*compuate_NoSlipBC_flux
 * functor to compute the contribution of a noSlip wall boundary condition
 * Sets the state for the flux evaluation such that the velocity at the
 * wall is zero.
 */
template<class InviscidFluxType, class ViscousFluxType>
struct compute_NoSlipBC_flux {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  cell_storage_field_type cell_flux_;
  vector_field_type cell_coordinates_;
  vector_field_type face_coordinates_, face_normal_, face_tangent_,
      face_binormal_;
  InviscidFluxType inviscid_flux_evaluator_;
  ViscousFluxType viscous_flux_evaluator_;

  compute_NoSlipBC_flux(Faces faces, solution_field_type cell_values,
      Cells cells, InviscidFluxType inviscid_flux,
      ViscousFluxType viscous_flux) :
          inviscid_flux_evaluator_(
          inviscid_flux), viscous_flux_evaluator_(viscous_flux) {

            std::copy(cell_values, cell_values + 5, cell_values_);
            std::copy(cells.cell_flux_, cells.cell_flux_ + 5, cell_flux_);

            for (int i = 0; i < 3; i++) {
              face_normal_[i]   = faces.face_normal_[i];
              face_tangent_[i]  = faces.face_tangent_[i];
              face_binormal_[i] = faces.face_binormal_[i];
              cell_coordinates_[i] = cells.coordinates_[i];
              face_coordinates_[i] = faces.coordinates_[i];
            }

            for (int i = 0; i < 2; i++) {
              face_cell_conn_[i] = faces.face_cell_conn_[i];
              cell_flux_index_[i] = faces.cell_flux_index_[i];
            }
  }

  void operator()(int i) const {
    int index = face_cell_conn_[0][i];

    double iflux[5];
    double vflux[5];
    double conservatives[5];
    double primitives_r[5];
    double primitives_l[5];

    for (int icomp = 0; icomp < 5; icomp++) {
      conservatives[icomp] = cell_values_[icomp][index];
      vflux[icomp] = 0.0;
    }

    ComputePrimitives(conservatives, primitives_l);

    //scale normal since it includes area.
    double area_norm = 0;
    for (int icomp = 0; icomp < 3; icomp++) {
        area_norm += face_normal_[icomp][i] * face_normal_[icomp][i];
    }
    area_norm = std::sqrt(area_norm);

    double uboundary = 0.0;
    uboundary += primitives_l[1] * face_normal_[0][i] / area_norm;
    uboundary += primitives_l[2] * face_normal_[1][i] / area_norm;
    uboundary += primitives_l[3] * face_normal_[2][i] / area_norm;

    primitives_r[0] = primitives_l[0];
    primitives_r[1] = primitives_l[1] - 2 * uboundary * face_normal_[0][i] / area_norm;
    primitives_r[2] = primitives_l[2] - 2 * uboundary * face_normal_[1][i] / area_norm;
    primitives_r[3] = primitives_l[3] - 2 * uboundary * face_normal_[2][i] / area_norm;
    primitives_r[4] = primitives_l[4];

    inviscid_flux_evaluator_.compute_flux(primitives_l, primitives_r, iflux,
        &face_normal_[0][i], &face_tangent_[0][i], &face_binormal_[0][i]);

    if (ViscousFluxType::isViscous) {
      double primitives_face[5];
      double gradients_face[5][3];
      double distance_to_wall = 0;
      double unit_normal[3];
      primitives_face[0] = primitives_l[0];
      primitives_face[1] = 0.0;
      primitives_face[2] = 0.0;
      primitives_face[3] = 0.0;
      primitives_face[4] = primitives_l[4];
      for (int idir = 0; idir < 3; idir++) {
        distance_to_wall += std::pow(
            face_coordinates_[idir][i] - cell_coordinates_[idir][index], 2);
        unit_normal[idir] = face_normal_[idir][i] / area_norm;
      }
      double inv_distance_to_wall = 1.0 / std::sqrt(distance_to_wall);

      for (int icomp = 0; icomp < 5; icomp++) {
        for (int idir = 0; idir < 3; idir++) {
          gradients_face[icomp][idir] = (primitives_face[icomp]
              - primitives_l[icomp]) * unit_normal[idir] * inv_distance_to_wall;
        }
      }
      viscous_flux_evaluator_.compute_flux(gradients_face, primitives_face,
          &face_normal_[0][i], vflux);
    }

#ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; ++icomp)
    {
      cell_flux_[index][cell_flux_index_[i][0]][icomp] = -iflux[icomp]+vflux[icomp];
    }
#endif
  }
};

#endif
