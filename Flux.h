#ifndef _INCLUDE_FLUX_H_
#define _INCLUDE_FLUX_H_

#include <cstdio>

#include "Faces.h"
#include "GasModel.h"
#include "MathToolsDevice.h"
#include "VanAlbadaLimiter.h"

/* compute_face_flux
 * functor to compute the internal face flux contributions.
 * Uses the templated Inviscid and Inviscid flux types to
 * compute the contribution. This functor organizes
 * the data to pass to the functions that compute the flux
 * and puts the flux contribution in the appropriate place
 * using either Gather-Sum or Atomics for thread safety.
 */
template <bool second_order, class InviscidFluxType, class ViscousFluxType>
struct compute_face_flux {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;
  typedef typename ViewTypes::gradient_field_type gradient_field_type;

  face_cell_conn_type face_cell_conn_;   // 2
  face_cell_conn_type cell_flux_index_;  // 2
  solution_field_type cell_values_;      // 5
  gradient_field_type cell_gradients_;   // 5, 3
  solution_field_type cell_limiters_;    // 5
  vector_field_type cell_coordinates_;   // 3
  cell_storage_field_type cell_flux_;    // 5
  vector_field_type face_coordinates_, face_normal_, face_tangent_,
      face_binormal_;  // 3
  const int* permute_vector_;
  InviscidFluxType inviscid_flux_evaluator_;
  ViscousFluxType viscous_flux_evaluator_;

  compute_face_flux(Faces faces, solution_field_type cell_values,
                    gradient_field_type cell_gradients,
                    solution_field_type cell_limiters, Cells cells,
                    InviscidFluxType inviscid_flux,
                    ViscousFluxType viscous_flux)
      : inviscid_flux_evaluator_(inviscid_flux),
        viscous_flux_evaluator_(viscous_flux),
        permute_vector_(faces.permute_vector_) {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2,
              face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2,
              cell_flux_index_);
    std::copy(cell_values, cell_values + 5, cell_values_);
    std::copy(cell_limiters, cell_limiters + 5, cell_limiters_);

    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++) cell_gradients_[i][j] = cell_gradients[i][j];
    }
  }

  void operator()(const int& ii, Faces faces, Cells cells) const {
    const int i = permute_vector_[ii];
    const int left_index = face_cell_conn_[0][i];
    const int right_index = face_cell_conn_[1][i];

    double flux[5];
    double conservatives_l[5];
    double conservatives_r[5];
    double primitives_l[5];
    double primitives_r[5];

    for (int icomp = 0; icomp < 5; icomp++) {
      conservatives_l[icomp] = cell_values_[icomp][left_index];
      conservatives_r[icomp] = cell_values_[icomp][right_index];
    }

    ComputePrimitives(conservatives_l, primitives_l);
    ComputePrimitives(conservatives_r, primitives_r);

    if (second_order) {
      // Extrapolation
      for (int icomp = 0; icomp < 5; icomp++) {
        double gradient_primitive_l_tmp = 0;
        double gradient_primitive_r_tmp = 0;

        for (int idir = 0; idir < 3; idir++) {
          gradient_primitive_l_tmp += (faces.coordinates_[idir][i] -
                                       cells.coordinates_[idir][left_index]) *
                                      cell_gradients_[icomp][idir][left_index];

          gradient_primitive_r_tmp += (faces.coordinates_[idir][i] -
                                       cells.coordinates_[idir][right_index]) *
                                      cell_gradients_[icomp][idir][right_index];
        }

        primitives_l[icomp] +=
            gradient_primitive_l_tmp * cell_limiters_[icomp][left_index];
        primitives_r[icomp] +=
            gradient_primitive_r_tmp * cell_limiters_[icomp][right_index];
      }

    }  // End of second order

    inviscid_flux_evaluator_.compute_flux(
        primitives_l, primitives_r, flux, &faces.face_normal_[0][i],
        &faces.face_tangent_[0][i], &faces.face_binormal_[0][i]);

    if (ViscousFluxType::isViscous) {
      double primitives_face[5];
      double gradients_face[5][3];

      for (int icomp = 0; icomp < 5; icomp++) {
        primitives_face[icomp] =
            0.5 * (primitives_l[icomp] + primitives_r[icomp]);

        for (int idir = 0; idir < 3; idir++) {
          gradients_face[icomp][idir] =
              0.5 * (cell_gradients_[icomp][idir][left_index] +
                     cell_gradients_[icomp][idir][right_index]);
        }
      }

      double vflux[5];
      viscous_flux_evaluator_.compute_flux(gradients_face, primitives_face,
                                           &faces.face_normal_[0][i], vflux);

      for (int icomp = 0; icomp < 5; icomp++) {
        flux[icomp] -= vflux[icomp];
      }
    }

#ifdef ATOMICS_FLUX
    for (int icomp = 0; icomp < 5; icomp++) {
      /*double* left_cell = &cell_flux_[left_index][0][icomp];
      Kokkos::atomic_add(left_cell, -flux[icomp]);
      double* right_cell = &cell_flux_[right_index][0][icomp];
      Kokkos::atomic_add(right_cell, flux[icomp]);*/
      cells.cell_flux_[icomp][left_index][0] -= flux[icomp];
      cells.cell_flux_[icomp][right_index][0] += flux[icomp];
    }
#endif

#ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; ++icomp) {
      cell_flux_[left_index][cell_flux_index_[i][0]][icomp] = -flux[icomp];
      cell_flux_[right_index][cell_flux_index_[i][1]][icomp] = flux[icomp];
    }
#endif
  }
};

/* apply_cell_flux
 * functor add the flux contributions to the residual
 * uses either gather-sum or atomics for thread safety
 */
struct apply_cell_flux {
  typedef typename ViewTypes::scalar_field_type scalar_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  typedef typename ViewTypes::cell_face_conn_type cell_face_conn_type;

  // int number_faces_;
  scalar_field_type volume_;
  // cell_storage_field_type flux_;
  solution_field_type residuals_;

  double dt_;
  double vol_;

  apply_cell_flux(double dt) : dt_(dt) {}

  void operator()(int i, Cells cells, solution_field_type res_vec) const {
    for (int icomp = 0; icomp < 5; icomp++) {
      res_vec[icomp][i] = 0.0;
    }

#ifdef ATOMICS_FLUX
    for (int flux_id = 0; flux_id < 1; flux_id++)
#else
    for (int flux_id = 0; flux_id < cells.nfaces_; flux_id++)
#endif
    {
      for (int icomp = 0; icomp < 5; icomp++) {
        res_vec[icomp][i] =
            res_vec[icomp][i] +
            dt_ / cells.volumes_[i] * cells.cell_flux_[icomp][i][flux_id];
      }
    }
  }
};

#endif
