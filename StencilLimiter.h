#include <algorithm>
#include <cstdio>

#include "Faces.h"
#include "MathToolsDevice.h"
#include "VenkatLimiter.h"
#include "ViewTypes.h"

#ifndef STENCILLIMITER_H_
#define STENCILLIMITER_H_

#define STENCIL_MAX(a, b)   \
  ({                        \
    __typeof__(a) _a = (a); \
    __typeof__(b) _b = (b); \
    _a > _b ? _a : _b;      \
  })
#define STENCIL_MIN(a, b)   \
  ({                        \
    __typeof__(a) _a = (a); \
    __typeof__(b) _b = (b); \
    _a < _b ? _a : _b;      \
  })

/* min_max_face
 * functor to compute the minimum and maximum value at each face
 * and scatters to the 2 connected elements.
 */
template <bool interior>
struct min_max_face {
  typedef typename ViewTypes::c_rnd_scalar_field_type scalar_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::c_rnd_vector_field_type vector_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  scalar_field_type cell_volumes_;
  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  vector_field_type face_normal_;
  cell_storage_field_type stencil_min_, stencil_max_;
  const int *permute_vector_;

  min_max_face(Faces faces, solution_field_type cell_values, Cells cells,
               cell_storage_field_type stencil_min,
               cell_storage_field_type stencil_max)
      :  // face_cell_conn_(faces.face_cell_conn_),
         // cell_flux_index_(faces.cell_flux_index_),
         // cell_values_(cell_values)
         // stencil_min_(stencil_min),
         // stencil_max_(stencil_max),
        permute_vector_(faces.permute_vector_) {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2,
              face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2,
              cell_flux_index_);
    std::copy(cell_values, cell_values + 5, cell_values_);
    std::copy(stencil_min, stencil_min + 5, stencil_min_);
    std::copy(stencil_max, stencil_max + 5, stencil_max_);
  }

  void operator()(const int &ii) const {
    const int i = permute_vector_[ii];

    const int left_index = face_cell_conn_[i][0];
    const int right_index = face_cell_conn_[i][1];

    double primitives_l[5];
    double primitives_r[5];

    const double gamma = 1.4;
    const double Rgas = 287.05;

    if (interior) {
      double r = cell_values_[0][left_index];
      double ri = 1.0 / r;
      double u = cell_values_[1][left_index] * ri;
      double v = cell_values_[2][left_index] * ri;
      double w = cell_values_[3][left_index] * ri;
      double k = 0.5 * (u * u + v * v + w * w);
      double e = cell_values_[4][left_index] * ri - k;
      double T = e * (gamma - 1.0) / Rgas;

      primitives_l[0] = r;
      primitives_l[1] = u;
      primitives_l[2] = v;
      primitives_l[3] = w;
      primitives_l[4] = T;

      r = cell_values_[0][right_index];
      ri = 1.0 / r;
      u = cell_values_[1][right_index] * ri;
      v = cell_values_[2][right_index] * ri;
      w = cell_values_[3][right_index] * ri;
      k = 0.5 * (u * u + v * v + w * w);
      e = cell_values_[4][right_index] * ri - k;
      T = e * (gamma - 1.0) / Rgas;

      primitives_r[0] = r;
      primitives_r[1] = u;
      primitives_r[2] = v;
      primitives_r[3] = w;
      primitives_r[4] = T;
    } else {
      const double r = cell_values_[0][left_index];
      const double ri = 1.0 / r;
      const double u = cell_values_[1][left_index] * ri;
      const double v = cell_values_[2][left_index] * ri;
      const double w = cell_values_[3][left_index] * ri;
      const double k = 0.5 * (u * u + v * v + w * w);
      const double e = cell_values_[4][left_index] * ri - k;
      const double T = e * (gamma - 1.0) / Rgas;

      primitives_l[0] = r;
      primitives_l[1] = u;
      primitives_l[2] = v;
      primitives_l[3] = w;
      primitives_l[4] = T;
    }

    const int cell_ind_0 = cell_flux_index_[0][i];
    const int cell_ind_1 = cell_flux_index_[1][i];

    for (int icomp = 0; icomp < 5; ++icomp) {
      const double face_min =
          interior ? STENCIL_MIN(primitives_r[icomp], primitives_l[icomp])
                   : primitives_l[icomp];
      const double face_max =
          interior ? STENCIL_MAX(primitives_r[icomp], primitives_l[icomp])
                   : primitives_l[icomp];

      /*#ifdef ATOMICS_FLUX
          //Need compare and exhange here instead of atomic add

          double * left_cell_min = &stencil_min_(left_index,0,icomp);
          bool success=false;
          do{
            double old_left_min =  *left_cell_min;
            double new_left_min = MathTools::min(*left_cell_min, face_min);
            double new_value =
      std::atomic_compare_exchange<double>(left_cell_min, old_left_min,
      new_left_min); success = new_value == new_left_min; } while(!success);
          double * left_cell_max = &stencil_max_(left_index,0,icomp);
          success=false;
          do{
            double old_left_max =  *left_cell_max;
            double new_left_max = MathTools::max(*left_cell_max, face_max);
            double new_value =
      std::atomic_compare_exchange<double>(left_cell_max, old_left_max,
      new_left_max); success = new_value == new_left_max; } while(!success);

          if(interior){
            double * right_cell_min = &stencil_min_(right_index,0,icomp);
            success=false;
            do{
              double old_right_min =  *right_cell_min;
              double new_right_min = MathTools::min(*right_cell_min, face_min);
              double new_value = Kokkos::atomic_compare_exchange(right_cell_min,
      old_right_min, new_right_min); success = new_value == new_right_min; }
      while(!success); double * right_cell_max =
      &stencil_max_(right_index,0,icomp); success=false; do{ double
      old_right_max =  *right_cell_max; double new_right_max =
      MathTools::max(*right_cell_max, face_max); double new_value =
      Kokkos::atomic_compare_exchange(right_cell_max, old_right_max,
      new_right_max); success = new_value == new_right_max; } while(!success);
          }
      #endif*/

#ifdef CELL_FLUX
      stencil_min_(left_index, cell_ind_0, icomp) = face_min;
      stencil_max_(left_index, cell_ind_0, icomp) = face_max;

      if (interior) {
        stencil_min_(right_index, cell_ind_1, icomp) = face_min;
        stencil_max_(right_index, cell_ind_1, icomp) = face_max;
      }
#endif
    }
  }
};

/* initialize_min_max
 * functor that initializes the minimum and maximum values using very large or
 * very small numbers.
 */
struct initialize_min_max {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  const int nfaces_;
  solution_field_type stencil_min_, stencil_max_;
  cell_storage_field_type stored_min_, stored_max_;

  initialize_min_max(int nfaces, solution_field_type stencil_min,
                     solution_field_type stencil_max,
                     cell_storage_field_type stored_min,
                     cell_storage_field_type stored_max)
      : nfaces_(nfaces) {
    std::copy(stencil_min, stencil_min + 5, stencil_min_);
    std::copy(stencil_max, stencil_max + 5, stencil_max_);
    std::copy(stored_min, stored_min + 5, stored_min_);
    std::copy(stored_max, stored_max + 5, stored_max_);
  }

  void operator()(int i) const {
    for (int icomp = 0; icomp < 5; ++icomp) {
      stencil_min_[icomp][i] = 1.0e300;
      stencil_max_[icomp][i] = -1.0e300;
      for (int iface = 0; iface < nfaces_; ++iface) {
        stored_min_[icomp][i][iface] = 1.0e300;
        stored_max_[icomp][i][iface] = -1.0e300;
      }
    }
  }
};

/* gather_min_max
 * functor that computes the minimum or maximum value of each variable over the
 * stencil of each cell.  The stencil consists of all cells which share a face
 * with this cell.
 */
struct gather_min_max {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  const int ncells_;
  const int nfaces_;
  cell_storage_field_type stored_min_, stored_max_;
  solution_field_type stencil_min_, stencil_max_;

  gather_min_max(Cells cells, cell_storage_field_type stored_min,
                 cell_storage_field_type stored_max,
                 solution_field_type stencil_min,
                 solution_field_type stencil_max)
      : ncells_(cells.ncells_), nfaces_(cells.nfaces_) {
    std::copy(stencil_min, stencil_min + 5, stencil_min_);
    std::copy(stencil_max, stencil_max + 5, stencil_max_);
    std::copy(stored_min, stored_min + 5, stored_min_);
    std::copy(stored_max, stored_max + 5, stored_max_);
  }
  void operator()(int i) const {
    for (int icomp = 0; icomp < 5; icomp++) {
#ifdef ATOMICS_FLUX
      stencil_min_[icomp][i] = stored_min_[icomp][i][0];
      stencil_max_[icomp][i] = stored_max_[icomp][i][0];
#endif

#ifdef CELL_FLUX
      for (int iface = 0; iface < nfaces_; ++iface) {
        stencil_min_[i][icomp] = MathTools::min(stencil_min_[i][icomp],
                                                stored_min_[i][iface][icomp]);
        stencil_max_[i][icomp] = MathTools::max(stencil_max_[i][icomp],
                                                stored_max_[i][iface][icomp]);
      }
#endif
    }
  }
};

/* initialize_limiter
 * functor that initializes the cell limiter value to 1.0.
 */
struct initialize_limiter {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  const int nfaces_;
  solution_field_type limiter_;
  cell_storage_field_type stored_limiter_;

  initialize_limiter(int nfaces, cell_storage_field_type stored_limiter,
                     solution_field_type limiter)
      : nfaces_(nfaces)
  // stored_limiter_(stored_limiter),
  // limiter_(limiter)
  {
    std::copy(stored_limiter, stored_limiter + 5, stored_limiter_);
    std::copy(limiter, limiter + 5, limiter_);
  }

  void operator()(int i) const {
    for (int icomp = 0; icomp < 5; icomp++) {
      for (int iface = 0; iface < nfaces_; iface++) {
        stored_limiter_[icomp][i][iface] = 1.0;
      }
      limiter_[icomp][i] = 1.0;
    }
  }
};

/* gather_limiter
 * functor that gathers and takes the minimum limiter value of the connected
 * faces.
 */
struct gather_limiter {
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;

  const int nfaces_;
  cell_storage_field_type stored_limiter_;
  solution_field_type limiter_;

  gather_limiter(int nfaces, cell_storage_field_type stored_limiter,
                 solution_field_type limiter)
      : nfaces_(nfaces) {
    std::copy(stored_limiter, stored_limiter + 5, stored_limiter_);
    std::copy(limiter, limiter + 5, limiter_);
  }

  void operator()(int i) const {
    for (int icomp = 0; icomp < 5; icomp++) {
#ifdef ATOMICS_FLUX
      limiter_[icomp][i] = stored_limiter_[icomp][i][0];
#endif

#ifdef CELL_FLUX
      for (int iface = 0; iface < nfaces_; ++iface) {
        limiter_[i][icomp] = MathTools::min(limiter_[i][icomp],
                                            stored_limiter_[i][iface][icomp]);
      }
#endif
    }
  }
};

/* limiter_face
 * functor that computes the limiter value for each face and scatter
 * contribution to the connected elements.  Uses gather-sum or atomics for
 * thread safety.
 */
template <bool interior>
struct limiter_face {
  typedef typename ViewTypes::c_rnd_scalar_field_type scalar_field_type;
  // typedef typename ViewTypes::c_rnd_solution_field_type solution_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  // typedef typename ViewTypes::c_rnd_face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::c_rnd_vector_field_type vector_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  // typedef typename ViewTypes::c_rnd_gradient_field_type gradient_field_type;
  typedef typename ViewTypes::gradient_field_type gradient_field_type;

  scalar_field_type cell_volumes_;
  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_min_, cell_max_, cell_values_;
  vector_field_type face_coordinates_, cell_coordinates_;
  gradient_field_type cell_gradients_;
  cell_storage_field_type limiter_;
  const int *permute_vector_;

  limiter_face(Faces faces, solution_field_type cell_values, Cells cells,
               gradient_field_type gradients, solution_field_type cell_min,
               solution_field_type cell_max, cell_storage_field_type limiter)
      :  // face_cell_conn_(faces.face_cell_conn_),
         // cell_flux_index_(faces.cell_flux_index_),
         // ell_min_(cell_min),
         // cell_max_(cell_max)
         // cell_values_(cell_values)
         // face_coordinates_(faces.coordinates_),
         // cell_coordinates_(cells.coordinates_)
         // cell_gradients_(gradients)
         // limiter_(limiter)
        permute_vector_(faces.permute_vector_) {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2,
              face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2,
              cell_flux_index_);
    std::copy(limiter, limiter + 5, limiter_);
    std::copy(cell_values, cell_values + 5, cell_values_);
    std::copy(cell_min, cell_min + 5, cell_min_);
    std::copy(cell_max, cell_max + 5, cell_max_);
    std::copy(faces.coordinates_, faces.coordinates_ + 5, face_coordinates_);
    std::copy(cells.coordinates_, cells.coordinates_ + 5, cell_coordinates_);

    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++) cell_gradients_[i][j] = gradients[i][j];
    }
  }

  void operator()(const int &ii) const {
    const int i = ii;  // permute_vector_[ii];
    const int left_index = face_cell_conn_[0][i];
    const int right_index = face_cell_conn_[1][i];

    double conservatives_l[5];
    double conservatives_r[5];
    double primitives_l[5];
    double primitives_r[5];

    for (int icomp = 0; icomp < 5; icomp++) {
      if (interior) {
        conservatives_l[icomp] = cell_values_[icomp][left_index];
        conservatives_r[icomp] = cell_values_[icomp][right_index];
      } else {
        conservatives_l[icomp] = cell_values_[icomp][left_index];
      }
    }

    if (interior) {
      ComputePrimitives(conservatives_l, primitives_l);
      ComputePrimitives(conservatives_r, primitives_r);
    } else {
      ComputePrimitives(conservatives_l, primitives_l);
    }

    // Compute left limiter value and compute right limiter value

    double limiter_left[5], limiter_right[5];
    // compute displacement and distance from cell center to face center.
    double displacement_l[3];
    double displacement_r[3];
    double distance_l = 0;
    double distance_r = 0;

    for (int idir = 0; idir < 3; idir++) {
      displacement_l[idir] =
          face_coordinates_[i][idir] - cell_coordinates_[left_index][idir];
      distance_l += displacement_l[idir] * displacement_l[idir];
      if (interior) {
        displacement_r[idir] =
            face_coordinates_[i][idir] - cell_coordinates_[right_index][idir];
        distance_r += displacement_r[idir] * displacement_r[idir];
      }
    }

    double dU_l[5];
    double dU_r[5];
    // Extrapolation
    for (int icomp = 0; icomp < 5; ++icomp) {
      dU_l[icomp] = 0;
      dU_r[icomp] = 0;
      for (int idir = 0; idir < 3; ++idir) {
        dU_l[icomp] +=
            displacement_l[idir] * cell_gradients_[left_index][icomp][idir];
        if (interior)
          dU_r[icomp] +=
              displacement_r[idir] * cell_gradients_[right_index][icomp][idir];
      }
    }

    for (int icomp = 0; icomp < 5; icomp++) {
      double dumax_l = cell_max_[left_index][icomp] - primitives_l[icomp];
      double dumin_l = cell_min_[left_index][icomp] - primitives_l[icomp];

      limiter_left[icomp] =
          VenkatLimiter::limit(dumax_l, dumin_l, dU_l[icomp], distance_l);
      if (interior) {
        double dumax_r = cell_max_[right_index][icomp] - primitives_r[icomp];
        double dumin_r = cell_min_[right_index][icomp] - primitives_r[icomp];
        limiter_right[icomp] =
            VenkatLimiter::limit(dumax_r, dumin_r, dU_r[icomp], distance_r);
      }
    }

    // Then write to memory
    /*#ifdef ATOMICS_FLUX
      for (int icomp = 0; icomp < 5; ++icomp)
      {
        double * left_cell_limiter = &limiter_[left_index][0][icomp];
        bool success=false;
        do{
          double old_left_limiter =  *left_cell_limiter;
          double new_left_limiter = MathTools::min(*left_cell_limiter,
    limiter_left[icomp]); double new_value =
    Kokkos::atomic_compare_exchange(left_cell_limiter, old_left_limiter,
    new_left_limiter); success = new_value == new_left_limiter; }
    while(!success);

        if(interior){
          double * right_cell_limiter = &limiter_(right_index,0,icomp);
          success=false;
          do{
            double old_right_limiter =  *right_cell_limiter;
            double new_right_limiter = MathTools::min(*right_cell_limiter,
    limiter_right[icomp]); double new_value =
    Kokkos::atomic_compare_exchange(right_cell_limiter, old_right_limiter,
    new_right_limiter); success = new_value == new_right_limiter; }
    while(!success);
        }
      }
    #endif*/

#ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; ++icomp) {
      limiter_[left_index][cell_flux_index_[i][0]][icomp] = limiter_left[icomp];

      if (interior) {
        limiter_[right_index][cell_flux_index_[i][1]][icomp] =
            limiter_right[icomp];
      }
    }
#endif
  }
};

/*StencilLimiter
 * Class to compute cell value of the stencil limiter.
 * Has methods to compute the min/max value of stencil
 * and compute the limiter value for each cell.
 */
class StencilLimiter {
  typedef typename ViewTypes::scalar_field_type scalar_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::cell_storage_field_type cell_storage_field_type;
  typedef typename ViewTypes::gradient_field_type gradient_field_type;

 public:
  StencilLimiter() {}
  StencilLimiter(Faces *internal_faces, std::vector<Faces *> *bc_faces,
                 Cells *cells, struct MeshData *mesh_data, int total_send_count,
                 int total_recv_count)
      : internal_faces_(internal_faces),
        bc_faces_(bc_faces),
        cells_(cells),
        mesh_data_(mesh_data)
  // ghosted_vars("ghosted_vars", total_recv_count*5),
  // ghosted_vars_host(ghosted_vars),
  // shared_vars("shared_vars", total_send_count*5),
  // shared_vars_host(shared_vars),
  // stored_min_("stored_min", cells->ncells_*5, cells->nfaces_),
  // stored_max_("stored_max", cells->ncells_*5, cells->nfaces_),
  // stored_limiter_("stored_limiter", cells->ncells_*5, cells->nfaces_)
  // stencil_min_("stencil_min", cells->ncells_*5),
  // stencil_max_("stencil_max", cells->ncells_*5)
  {
    for (int i = 0; i < 5; i++) {
      stencil_min_[i] = (double *)malloc(sizeof(double) * cells->ncells_ * 5);
      stencil_max_[i] = (double *)malloc(sizeof(double) * cells->ncells_ * 5);

      stored_min_[i] = (double **)new double[cells->ncells_ * 5];
      stored_max_[i] = (double **)new double[cells->ncells_ * 5];
      stored_limiter_[i] = new double *[cells->ncells_];

      for (int j = 0; j < cells->ncells_; j++) {
        // stored_min_[i][j] = new double[cells->nfaces_];
        // stored_max_[i][j] = new double[cells->nfaces_];
        stored_limiter_[i][j] = new double[cells->nfaces_];
      }
    }
  }

  void compute_min_max(solution_field_type sol_np1_vec) {
    initialize_min_max init_min_max(cells_->nfaces_, stencil_min_, stencil_max_,
                                    stored_min_, stored_max_);
    // Kokkos::parallel_for(mesh_data_->num_owned_cells, init_min_max);

    // Internal Faces
    const int ninternal_faces = internal_faces_->nfaces_;
    min_max_face<true> min_max_internal(*internal_faces_, sol_np1_vec, *cells_,
                                        stored_min_, stored_max_);
    // Kokkos::parallel_for(ninternal_faces, min_max_internal);

    // Boundary Faces
    typename std::vector<Faces *>::iterator bcf_iter, bcf_iter_end;
    bcf_iter = bc_faces_->begin();
    bcf_iter_end = bc_faces_->end();

    for (; bcf_iter != bcf_iter_end; ++bcf_iter) {
      Faces *faces = *bcf_iter;
      const int nboundary_faces = faces->nfaces_;
      min_max_face<false> bc_min_max(*faces, sol_np1_vec, *cells_, stored_min_,
                                     stored_max_);
      // Kokkos::parallel_for(nboundary_faces, bc_min_max);
    }

    // Kokkos::fence();

    gather_min_max gather(*cells_, stored_min_, stored_max_, stencil_min_,
                          stencil_max_);
    // Kokkos::parallel_for(mesh_data_->num_owned_cells, gather);
    // Kokkos::fence();
  }

  void communicate_min_max() {
    /*
  // For min
      extract_shared_vector<5> extract_shared_min(stencil_min_,
  mesh_data_->send_local_ids, shared_vars);
      // Kokkos::parallel_for(mesh_data_->num_ghosts, extract_shared_min);
      // Kokkos::fence();
      // Kokkos::deep_copy(shared_vars_host, shared_vars);

      communicate_ghosted_cell_data(mesh_data_->sendCount,
  mesh_data_->recvCount,
  shared_vars_host.ptr_on_device(),ghosted_vars_host.ptr_on_device(), 5);

      // Kokkos::deep_copy(ghosted_vars, ghosted_vars_host);
      insert_ghost_vector<5> insert_ghost_min(stencil_min_,
  mesh_data_->recv_local_ids, ghosted_vars);
      // Kokkos::parallel_for(mesh_data_->num_ghosts, insert_ghost_min);
      // Kokkos::fence();

  // For max
      extract_shared_vector<5> extract_shared_max(stencil_max_,
  mesh_data_->send_local_ids, shared_vars);
      // Kokkos::parallel_for(mesh_data_->num_ghosts, extract_shared_max);
      // Kokkos::fence();
      // Kokkos::deep_copy(shared_vars_host, shared_vars);

      communicate_ghosted_cell_data(mesh_data_->sendCount,
  mesh_data_->recvCount,
  shared_vars_host.ptr_on_device(),ghosted_vars_host.ptr_on_device(), 5);

      // Kokkos::deep_copy(ghosted_vars, ghosted_vars_host);
      insert_ghost_vector<5> insert_ghost_max(stencil_max_,
  mesh_data_->recv_local_ids, ghosted_vars);
      // Kokkos::parallel_for(mesh_data_->num_ghosts, insert_ghost_max);
      // Kokkos::fence();
  // TODO: Maybe combined or overlapped in future.*/
  }

  void compute_limiter(solution_field_type sol_np1_vec,
                       solution_field_type limiter,
                       gradient_field_type gradients) {
    initialize_limiter init_limiter(cells_->nfaces_, stored_limiter_, limiter);
    for (int i = 0; i < mesh_data_->num_owned_cells; i++) {
      init_limiter(i);
    }
    // Kokkos::parallel_for(mesh_data_->num_owned_cells, init_limiter);

    // Internal Faces
    const int ninternal_faces = internal_faces_->nfaces_;
    limiter_face<true> limiter_internal(*internal_faces_, sol_np1_vec, *cells_,
                                        gradients, stencil_min_, stencil_max_,
                                        stored_limiter_);
    for (int i = 0; i < ninternal_faces; i++) limiter_internal(i);
    // Kokkos::parallel_for(ninternal_faces, limiter_internal);

    // Boundary Faces
    typename std::vector<Faces *>::iterator bcf_iter, bcf_iter_end;
    bcf_iter = bc_faces_->begin();
    bcf_iter_end = bc_faces_->end();
    for (; bcf_iter != bcf_iter_end; ++bcf_iter) {
      Faces *faces = *bcf_iter;
      const int nboundary_faces = faces->nfaces_;
      limiter_face<false> limiter_bc(*faces, sol_np1_vec, *cells_, gradients,
                                     stencil_min_, stencil_max_,
                                     stored_limiter_);
      for (int i = 0; i < nboundary_faces; i++) limiter_bc(i);
      // Kokkos::parallel_for(nboundary_faces, limiter_bc);
    }
    // Kokkos::fence();

    gather_limiter gather(cells_->nfaces_, stored_limiter_, limiter);
    for (int i = 0; i < mesh_data_->num_owned_cells; i++) gather(i);
    // Kokkos::parallel_for(mesh_data_->num_owned_cells, gather);
    // Kokkos::fence();
  }

  /*void communicate_limiter(solution_field_type limiter) {

      extract_shared_vector<5> extract_shared_limiter(limiter,
  mesh_data_->send_local_ids, shared_vars);
      // Kokkos::parallel_for(mesh_data_->num_ghosts, extract_shared_limiter);
      // Kokkos::fence();
      // Kokkos::deep_copy(shared_vars_host, shared_vars);

      communicate_ghosted_cell_data(mesh_data_->sendCount,
  mesh_data_->recvCount, shared_vars_host.ptr_on_device(),
  ghosted_vars_host.ptr_on_device(), 5);

      // Kokkos::deep_copy(ghosted_vars, ghosted_vars_host);
      insert_ghost_vector<5> insert_ghost_limiter(limiter,
  mesh_data_->recv_local_ids, ghosted_vars);
      // Kokkos::parallel_for(mesh_data_->num_ghosts, insert_ghost_limiter);
      // Kokkos::fence();
  }*/

 private:
  Faces *internal_faces_;
  std::vector<Faces *> *bc_faces_;
  Cells *cells_;
  struct MeshData *mesh_data_;
  // scalar_field_type ghosted_vars;
  // typename scalar_field_type::HostMirror ghosted_vars_host;
  // scalar_field_type shared_vars;
  // typename scalar_field_type::HostMirror shared_vars_host;
  cell_storage_field_type stored_min_;
  cell_storage_field_type stored_max_;
  cell_storage_field_type stored_limiter_;
  solution_field_type stencil_min_;
  solution_field_type stencil_max_;
};

#endif
