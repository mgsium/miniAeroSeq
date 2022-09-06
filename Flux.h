#ifndef _INCLUDE_FLUX_H_
#define _INCLUDE_FLUX_H_

#include <cstdio>
#include "MathToolsDevice.h"
#include "Faces.h"
#include "VanAlbadaLimiter.h"
#include "GasModel.h"


/* compute_face_flux
 * functor to compute the internal face flux contributions.
 * Uses the templated Inviscid and Inviscid flux types to
 * compute the contribution. This functor organizes
 * the data to pass to the functions that compute the flux
 * and puts the flux contribution in the appropriate place
 * using either Gather-Sum or Atomics for thread safety.
 */
template<bool second_order, class InviscidFluxType, class ViscousFluxType >
struct compute_face_flux {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;
  typedef typename ViewTypes::gradient_field_type gradient_field_type;

  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  gradient_field_type cell_gradients_;
  solution_field_type cell_limiters_;
  vector_field_type cell_coordinates_;
  cell_storage_field_type cell_flux_;
  vector_field_type face_coordinates_, face_normal_, face_tangent_,
      face_binormal_;
  const int* permute_vector_;
  InviscidFluxType inviscid_flux_evaluator_;
  ViscousFluxType viscous_flux_evaluator_;

  compute_face_flux(Faces faces, solution_field_type cell_values,
      gradient_field_type cell_gradients, solution_field_type cell_limiters,
      Cells cells, InviscidFluxType inviscid_flux,
      ViscousFluxType viscous_flux) :
      /*face_cell_conn_(faces.face_cell_conn_), cell_flux_index_(
          faces.cell_flux_index_), cell_values_(cell_values),*/ 
          // cell_gradients_(cell_gradients), 
          // cell_limiters_(cell_limiters), 
          // cell_coordinates_(cells.coordinates_), 
          // cell_flux_(cells.cell_flux_), 
          // face_coordinates_(faces.coordinates_), 
          // face_normal_(faces.face_normal_), 
          // face_tangent_(faces.face_tangent_), 
          // face_binormal_(faces.face_binormal_), 
          inviscid_flux_evaluator_(inviscid_flux), 
          viscous_flux_evaluator_(viscous_flux) , 
          permute_vector_(faces.permute_vector_) {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2, face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2, cell_flux_index_);
    std::copy(cell_values, cell_values + 5, cell_values);
    std::copy(cell_limiters, cell_limiters + 5, cell_limiters_);
    std::copy(cells.coordinates_, cells.coordinates_ + 3, cell_coordinates_);
    std::copy(cells.cell_flux_, cells.cell_flux_ + 5, cell_flux_);
    std::copy(faces.coordinates_, faces.coordinates_ + 3, face_coordinates_);
    std::copy(faces.face_normal_, faces.face_normal_ + 3, face_normal_);
    std::copy(faces.face_tangent_, faces.face_tangent_ + 3, face_tangent_);
    std::copy(faces.face_binormal_, faces.face_binormal_ + 3, face_binormal_);

    for(int i = 0; i < 5; i++) {
      for(int j = 0; j < 3; j++)
        cell_gradients_[i][j] = cell_gradients[i][j];
    }
  }
  
  void operator()(const int& ii) const {
    const int i = permute_vector_[ii];
    const int left_index = face_cell_conn_[i][0];
    const int right_index = face_cell_conn_[i][1];

    double flux[5];
    double conservatives_l[5];
    double conservatives_r[5];
    double primitives_l[5];
    double primitives_r[5];

    for (int icomp = 0; icomp < 5; ++icomp) {
      conservatives_l[icomp] = cell_values_[left_index][icomp];
      conservatives_r[icomp] = cell_values_[right_index][icomp];
    }

    ComputePrimitives(conservatives_l, primitives_l);
    ComputePrimitives(conservatives_r, primitives_r);

    if (second_order) {

      //Extrapolation
      for (int icomp = 0; icomp < 5; ++icomp) {
      	double gradient_primitive_l_tmp = 0;
	      double gradient_primitive_r_tmp = 0;

        for (int idir = 0; idir < 3; ++idir) {
    	    gradient_primitive_l_tmp += (face_coordinates_[i][idir]
            	- cell_coordinates_[left_index][idir])
		          * cell_gradients_[left_index][icomp][idir];

    	    gradient_primitive_r_tmp += (face_coordinates_[i][idir]
		          - cell_coordinates_[right_index][idir])
		          * cell_gradients_[right_index][icomp][idir];
        }

        primitives_l[icomp] += gradient_primitive_l_tmp *
                cell_limiters_[left_index][icomp];
        primitives_r[icomp] += gradient_primitive_r_tmp *
                cell_limiters_[right_index][icomp];
      }

    } // End of second order


    inviscid_flux_evaluator_.compute_flux(primitives_l, primitives_r, flux,
        &face_normal_[i][0], &face_tangent_[i][0], &face_binormal_[i][0]);

    if (ViscousFluxType::isViscous) {
      double primitives_face[5];
      double gradients_face[5][3];

      for (int icomp = 0; icomp < 5; ++icomp) {
        primitives_face[icomp] = 0.5
            * (primitives_l[icomp] + primitives_r[icomp]);

        for (int idir = 0; idir < 3; ++idir) {
          gradients_face[icomp][idir] = 0.5
              * (cell_gradients_[left_index][icomp][idir]
                  + cell_gradients_[right_index][icomp][idir]);
        }
      }

      double vflux[5];
      viscous_flux_evaluator_.compute_flux(gradients_face, primitives_face,
          &face_normal_[i][0], vflux);

      for (int icomp = 0; icomp < 5; ++icomp) {
        flux[icomp] -= vflux[icomp];
      }
    }

#ifdef ATOMICS_FLUX
    for (int icomp = 0; icomp < 5; ++icomp)
    {
      double * left_cell = &cell_flux_(left_index,0,icomp);
      Kokkos::atomic_add(left_cell, -flux[icomp]);
      double * right_cell = &cell_flux_(right_index,0,icomp);
      Kokkos::atomic_add(right_cell, flux[icomp]);
    }
#endif

#ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; ++icomp)
    {
      cell_flux_(left_index,cell_flux_index_(i,0),icomp) = -flux[icomp];
      cell_flux_(right_index,cell_flux_index_(i,1),icomp) = flux[icomp];
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

  int number_faces_;
  scalar_field_type volume_;
  cell_storage_field_type flux_;
  solution_field_type residuals_;

  double dt_;
  double vol_;

  apply_cell_flux(Cells cells, solution_field_type residuals, double dt) :
      number_faces_(cells.nfaces_), 
      volume_(cells.volumes_),
      // flux_(cells.cell_flux_), 
      // residuals_(residuals), 
      dt_(dt) {
    std::copy(cells.cell_flux_, cells.cell_flux_ + 5, flux_);
    std::copy(residuals, residuals + 5, residuals_);
  }

  void operator()(int i) const {

    for (int icomp = 0; icomp < 5; ++icomp) {
      residuals_[i][icomp] = 0.0;
    }
#ifdef ATOMICS_FLUX
    for(int flux_id=0; flux_id<1; ++flux_id)
#else
    for (int flux_id = 0; flux_id < number_faces_; ++flux_id)
#endif

        {
      for (int icomp = 0; icomp < 5; ++icomp) {
        residuals_[i][icomp] = residuals_[i][icomp]
            + dt_ / volume_[i] * flux_[i][flux_id][icomp];
      }
    }
  }
};

#endif