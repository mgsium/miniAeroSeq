#ifndef INCLUDE_CELLS_H_
#define INCLUDE_CELLS_H_
#include <memory>

#include "Cell.h"
#include "ViewTypes.h"

/*Cells
 * struct containing the cell data used for simulation.
 * Includes things such as coordinates, volume, the gradient, and
 * the flux contributions to the residual.
 */
struct Cells {
  typedef typename ViewTypes::scalar_field_type scalar_field_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::cell_face_conn_type cell_face_conn_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  typedef typename ViewTypes::gradient_storage_field_type
      gradient_storage_field_type;

 public:
  int ncells_;
  int nfaces_;
  vector_field_type coordinates_;
  scalar_field_type volumes_;
  cell_storage_field_type cell_flux_;
  gradient_storage_field_type cell_gradient_;

  Cells() {}

  Cells(int ncells, int faces_per_elem)
      : ncells_(ncells),
        nfaces_(faces_per_elem)
  // coordinates_("cell_coordinates", ncells)
  // volumes_("cell_volumes", ncells)
  // TODO: Make this more efficient for atomics.
  {
    // Allocating memory for volumes
    volumes_ = new double[ncells];

    int mult = faces_per_elem;
#ifndef ATOMICS_FLUX
    mult = faces_per_elem;
#endif

    // Allocating memory for flux & coordinates
    for (int i = 0; i < 5; i++) {
      cell_flux_[i] = new double*[ncells];
      for (int j = 0; j < ncells; j++) cell_flux_[i][j] = new double[mult];
    }

    for (int i = 0; i < 3; i++) coordinates_[i] = new double[ncells];

    // Allocating memory for gradient storage
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++) {
        cell_gradient_[i][j] = new double*[ncells];
        for (int z = 0; z < ncells; z++)
          cell_gradient_[i][j][z] = new double[mult];
      }
    }
  }
};

/*zero_cell_flux
 * Functor to reset the flux contributions to the residual
 * to zero.
 */
struct zero_cell_flux {
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  typedef typename ViewTypes::gradient_storage_field_type
      gradient_storage_field_type;

  const int ncells_;
  const int nfaces_;
  // Cells cells_;
  // cell_storage_field_type cell_flux_;
  gradient_storage_field_type cell_gradient_;

  zero_cell_flux(Cells cells)
      : ncells_(cells.ncells_), nfaces_(cells.nfaces_) {}

  void operator()(int i, Cells cells) const {
    for (int icomp = 0; icomp < 5; icomp++) {
#ifdef ATOMICS_FLUX
      cells.cell_flux_[icomp][i][0] = 0.0;
#else
      for (int iface = 0; iface < cells.nfaces_; iface++) {
        cells.cell_flux_[icomp][i][iface] = 0.0;
      }
#endif
      for (int iDir = 0; iDir < 3; iDir++) {
        cells.cell_gradient_[icomp][iDir][i][0] = 0.0;
      }
    }
  }
};

/*copy_cell_data
 * functor to copy the Cell information from the setup datastructure to
 * Kokkos datastructure.
 */
/*void copy_cell_data(Cells device_cells, std::vector<Cell> & mesh_cells){

  typedef double * scalar_field_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;

  int ncells = mesh_cells.size();

    for(int i = 0; i < mesh_cells.size(); ++i){
      device_cells.volumes_[i] = mesh_cells[i].GetVolume();
      for(int j=0; j<3; ++j){
        double * tmp_coord = mesh_cells[i].GetCoords();
        device_cells.coordinates_[i][j] = tmp_coord[j];
      }
    }
}*/

#endif