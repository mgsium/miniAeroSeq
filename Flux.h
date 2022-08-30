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
template<bool second_order, class InviscidFluxType, class ViscousFluxType>
void compute_face_flux(
    Faces faces,
    ViewTypes::solution_field_type cell_values_,
    ViewTypes::gradient_field_type cell_gradients_,
    ViewTypes::solution_field_type cell_limiters_,
    Cells cells, 
    InviscidFluxType inviscid_flux,
    ViscousFluxType viscous_flux
) {
    ViewTypes::face_cell_conn_type face_cell_conn_ = faces.face_cell_conn_;
    ViewTypes::face_cell_conn_type cell_flux_index_ = faces.cell_flux_index_;
    ViewTypes::vector_field_type cell_coordinates_ = cells.coordinates_;
    ViewTypes::cell_storage_field_type cell_flux_ = cells.cell_flux_;
    ViewTypes::vector_field_type face_coordinates_ = faces.coordinates_;
    ViewTypes::vector_field_type face_normal_ = faces.face_normal_;
    ViewTypes::vector_field_type face_tangent_ = faces.face_tangent_;
    ViewTypes::vector_field_type face_binormal_ = faces.face_binormal_;

    const int i = permute_vector_[ii];
    const int left_index = face_cell_conn_[i][0]];
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

    ComputePrimitives[conservatives_l][primitives_l];
    ComputePrimitives[conservatives_r][primitives_r];

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
}

#endif