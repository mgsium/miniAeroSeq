#ifndef _INCLUDE_EXTRAPOLATE_BC_H_
#define _INCLUDE_EXTRAPOLATE_BC_H_

#include <cstdio>
#include "Faces.h"
#include "GasModel.h"

/* compute_extrapolate_BC_flux
 * functor to compute the contribution of an extrapolation boundary condition
 * state is extrapolated using zeroth order extrapolation(copy) to external cell
 * and external face flux contribution is then computed.
 */
template<class FluxType>
void compute_extrapolateBC_flux (
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

};

#endif