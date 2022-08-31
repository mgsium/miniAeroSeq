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
      face_cell_conn_(faces.face_cell_conn_), cell_flux_index_(
          faces.cell_flux_index_), cell_values_(cell_values), cell_flux_(
          cells.cell_flux_), cell_coordinates_(cells.coordinates_), face_coordinates_(
          faces.coordinates_), face_normal_(faces.face_normal_), face_tangent_(
          faces.face_tangent_), face_binormal_(faces.face_binormal_), inviscid_flux_evaluator_(
          inviscid_flux), viscous_flux_evaluator_(viscous_flux) {
  }

  void operator()(int i) const {
    int index = face_cell_conn_[i][0];

    double iflux[5];
    double vflux[5];
    double conservatives[5];
    double primitives_r[5];
    double primitives_l[5];

    for (int icomp = 0; icomp < 5; ++icomp) {
      conservatives[icomp] = cell_values_[index][icomp];
      vflux[icomp] = 0.0;
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

    inviscid_flux_evaluator_.compute_flux(primitives_l, primitives_r, iflux,
        &face_normal_[i][0], &face_tangent_[i][0], &face_binormal_[i][0]);

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
      for (int idir = 0; idir < 3; ++idir) {
        distance_to_wall += std::pow(
            face_coordinates_[i][idir] - cell_coordinates_[index][idir], 2);
        unit_normal[idir] = face_normal_[i][idir] / area_norm;
      }
      double inv_distance_to_wall = 1.0 / std::sqrt(distance_to_wall);

      for (int icomp = 0; icomp < 5; ++icomp) {
        for (int idir = 0; idir < 3; ++idir) {
          gradients_face[icomp][idir] = (primitives_face[icomp]
              - primitives_l[icomp]) * unit_normal[idir] * inv_distance_to_wall;
        }
      }
      viscous_flux_evaluator_.compute_flux(gradients_face, primitives_face,
          &face_normal_[i][0], vflux);
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
