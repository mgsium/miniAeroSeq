#ifndef _INCLUDE_EXTRAPOLATE_BC_H_
#define _INCLUDE_EXTRAPOLATE_BC_H_

#include <algorithm>
#include <cstdio>
#include "Faces.h"
#include "GasModel.h"

/* compute_extrapolate_BC_flux
 * functor to compute the contribution of an extrapolation boundary condition
 * state is extrapolated using zeroth order extrapolation(copy) to external cell
 * and external face flux contribution is then computed.
 */
template<class FluxType>
struct compute_extrapolateBC_flux {
  typedef typename ViewTypes::c_rnd_solution_field_type solution_field_type;
  typedef typename ViewTypes::c_rnd_face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::c_vector_field_type vector_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  cell_storage_field_type cell_flux_;
  vector_field_type face_normal_, face_tangent_, face_binormal_;
  FluxType flux_evaluator_;

  compute_extrapolateBC_flux(Faces faces,
      solution_field_type cell_values, Cells cells, FluxType flux) :
      flux_evaluator_(flux) {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2, face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2, cell_flux_index_);
    std::copy(cell_values, cell_values + 5, cell_values);
    std::copy(cells.cell_flux_, cells.cell_flux_ + 5, cell_flux_);
    std::copy(faces.face_normal_, faces.face_normal_ + 3, face_normal_);
    std::copy(faces.face_tangent_, faces.face_tangent_ + 3, face_tangent_);
    std::copy(faces.face_binormal_, faces.face_binormal_ + 3, face_binormal_);
  }

  void operator()(int i) const {
    int index = face_cell_conn_[i][0];

    double flux[5];
    double conservatives[5];
    double primitives[5];

    for (int icomp = 0; icomp < 5; ++icomp) {
        conservatives[icomp] = cell_values_[index][icomp];
    }

    ComputePrimitives(conservatives, primitives);

    flux_evaluator_.compute_flux(primitives, primitives, flux, &face_normal_[i][0],
        &face_tangent_[i][0], &face_binormal_[i][0]);

    /*#ifdef ATOMICS_FLUX
    for (int icomp = 0; icomp < 5; ++icomp)
    {
        double * cell = &cell_flux_(index,0,icomp);
        Kokkos::atomic_add(cell, -flux[icomp]);
    }
    #endif*/

    #ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; ++icomp)
    {
        cell_flux_[index][cell_flux_index_[i][0]][icomp] = -flux[icomp];
    }
    #endif
  }
};

#endif