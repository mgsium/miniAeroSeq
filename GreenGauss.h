#include <cstdio>
#include "Faces.h"

#ifndef GREENGAUSS_H_
#define GREENGAUSS_H_

/*green_gauss_face
 * functor to compute internal face contributions for Green-Gauss gradient computation.
 */
struct green_gauss_face{
  typedef typename ViewTypes::c_rnd_scalar_field_type scalar_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::c_rnd_face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::c_rnd_vector_field_type vector_field_type;
  typedef typename ViewTypes::gradient_storage_field_type gradient_storage_field_type;

  scalar_field_type cell_volumes_;
  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  gradient_storage_field_type cell_gradient_;
  vector_field_type face_normal_;
  const int* permute_vector_;

  green_gauss_face(Faces faces, solution_field_type cell_values, Cells cells):
    cell_volumes_(cells.volumes_),
    // cell_gradient_(cells.cell_gradient_)
    permute_vector_(faces.permute_vector_)
  {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2, face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2, cell_flux_index_);
    std::copy(cell_values, cell_values + 5, cell_values);
    std::copy(faces.face_normal_, faces.face_normal_ + 3, face_normal_);

    for(int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++)
        cell_gradient_[i][j] = cells.cell_gradient_[i][j];
    }
  }

  void operator()( const int& ii )const{
    const int i = permute_vector_[ii];
    const int left_index = face_cell_conn_[i][0];
    const int right_index = face_cell_conn_[i][1];

    const double gamma = 1.4;
    const double Rgas = 287.05;

    const double left_r  = cell_values_[left_index][0];
    const double left_ri = 1.0 / left_r;
    const double left_u  = cell_values_[left_index][1] * left_ri;
    const double left_v  = cell_values_[left_index][2] * left_ri;
    const double left_w  = cell_values_[left_index][3] * left_ri;
    const double left_k  = 0.5 * (left_u * left_u + left_v * left_v + left_w * left_w);
    const double left_e  = cell_values_[left_index][4] * left_ri - left_k;
    const double left_T  = left_e * (gamma - 1.0) / Rgas;

    const double primitives_l[5] = { left_r, left_u, left_v, left_w, left_T };

    const double right_r  = cell_values_[right_index][0];
    const double right_ri = 1.0 / right_r;
    const double right_u  = cell_values_[right_index][1] * right_ri;
    const double right_v  = cell_values_[right_index][2] * right_ri;
    const double right_w  = cell_values_[right_index][3] * right_ri;
    const double right_k  = 0.5 * (right_u * right_u + right_v * right_v + right_w * right_w);
    const double right_e  = cell_values_[right_index][4] * right_ri - right_k;
    const double right_T  = right_e * (gamma - 1.0) / Rgas;

    const double primitives_r[5] = { right_r, right_u, right_v, right_w, right_T };

    const double cell_vol_left = cell_volumes_[left_index];
    const double cell_vol_right = cell_volumes_[right_index];
    const int cell_ind_0 = cell_flux_index_[i][0];
    const int cell_ind_1 = cell_flux_index_[i][1];

    for(int idir = 0; idir < 3; ++idir)
    {
        const double face_norm = face_normal_[i][idir];

        for(int icomp = 0; icomp < 5; ++icomp) {
            const double gradient = 0.5*(primitives_l[icomp]+primitives_r[icomp])*face_norm;

/*#ifdef ATOMICS_FLUX
            double * left_cell = &cell_gradient_(left_index,0,icomp,idir);
            Kokkos::atomic_add(left_cell, gradient/cell_vol_left);

            double * right_cell = &cell_gradient_(right_index,0,icomp,idir);
            Kokkos::atomic_add(right_cell, -gradient/cell_vol_right);
#endif*/

#ifdef CELL_FLUX
            cell_gradient_[left_index][cell_ind_0][icomp][idir] = gradient/cell_vol_left;
            cell_gradient_[right_index][cell_ind_1][icomp][idir] = -gradient/cell_vol_right;
#endif
        }
    }
  }
};

/*green_gauss_boundary_face
 * functor to compute boundary face contributions for Green-Gauss gradient computation.
 */
struct green_gauss_boundary_face{
  typedef typename ViewTypes::scalar_field_type scalar_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;
  typedef typename ViewTypes::gradient_storage_field_type gradient_storage_field_type;

  scalar_field_type cell_volumes_;
  face_cell_conn_type face_cell_conn_;
  face_cell_conn_type cell_flux_index_;
  solution_field_type cell_values_;
  gradient_storage_field_type cell_gradient_;
  vector_field_type face_normal_;

  green_gauss_boundary_face(Faces faces, solution_field_type cell_values, Cells cells):
    cell_volumes_(cells.volumes_)
    // cell_gradient_(cells.cell_gradient_)
  {
    std::copy(faces.face_cell_conn_, faces.face_cell_conn_ + 2, face_cell_conn_);
    std::copy(faces.cell_flux_index_, faces.cell_flux_index_ + 2, cell_flux_index_);
    std::copy(cell_values, cell_values + 5, cell_values);
    std::copy(faces.face_normal_, faces.face_normal_ + 3, face_normal_);

    for(int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++)
        cell_gradient_[i][j] = cells.cell_gradient_[i][j];
    }
  }

  void operator()( int i )const{
    int index = face_cell_conn_[i][0];

  const double gamma = 1.4;
  const double Rgas = 287.05;

  double gradient[5][3];
  double primitives[5];

    const double r  = cell_values_[index][0];
    const double ri = 1.0 / r;
    const double u  = cell_values_[index][1] * ri;
    const double v  = cell_values_[index][2] * ri;
    const double w  = cell_values_[index][3] * ri;
    const double k  = 0.5 * (u * u + v * v + w * w);
    const double e  = cell_values_[index][4] * ri - k;
    const double T  = e * (gamma - 1.0) / Rgas;

    primitives[0] = r;
    primitives[1] = u;
    primitives[2] = v;
    primitives[3] = w;
    primitives[4] = T;

  for(int icomp = 0; icomp < 5; ++icomp)
  {
    for(int idir = 0; idir < 3; ++idir)
    {
        gradient[icomp][idir] = primitives[icomp]*face_normal_[i][idir];
    }
  }

/*#ifdef ATOMICS_FLUX
  for (int icomp = 0; icomp < 5; ++icomp)
  {
    for(int idir = 0; idir < 3; ++idir)
    {
      double * cell = &cell_gradient_(index,0,icomp,idir);
      Kokkos::atomic_fetch_add(cell, gradient[icomp][idir]/cell_volumes_(index));
    }
  }
#endif*/

#ifdef CELL_FLUX
    for (int icomp = 0; icomp < 5; ++icomp)
    {
    for(int idir = 0; idir < 3; ++idir)
    {
        cell_gradient_[index][cell_flux_index_[i][0]][icomp][idir] = gradient[icomp][idir]/cell_volumes_[index];
    }
    }
#endif
}
};

/* green_gauss_gradient_sum
 * functor to sum all of the contributions to the gradient
 * uses either gather-sum or atomics for thread safety.
 */
struct green_gauss_gradient_sum{
  typedef typename ViewTypes::gradient_storage_field_type gradient_storage_field_type;
  typedef typename ViewTypes::gradient_field_type gradient_field_type;

  gradient_storage_field_type face_gradient_;
  gradient_field_type gradient_;
  int number_faces_;

  green_gauss_gradient_sum(Cells cells, gradient_field_type gradient):
    // face_gradient_(cells.cell_gradient_),
    // gradient_(gradient),
    number_faces_(cells.nfaces_)
  {
    for(int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++) {
        face_gradient_[i][j] = cells.cell_gradient_[i][j];
        gradient_[i][j] = gradient[i][j];
      }
    }
  }

  void operator()( int i )const{
    for(int icomp = 0; icomp < 5; ++icomp)
  {

    for(int idir = 0; idir < 3; ++idir)
    {
      gradient_[icomp][idir][i] = 0;
    }
  }
#ifdef ATOMICS_FLUX
    for(int iface=0; iface<1; ++iface)
#else
    for(int iface=0; iface<number_faces_; ++iface)
#endif
    {

      for(int icomp = 0; icomp < 5; ++icomp)
      {
        for(int idir = 0; idir < 3; ++idir)
        {
          gradient_[icomp][idir][i] += face_gradient_[iface][icomp][idir][i];
        }
      }
    }
  }
};

/*GreenGauss
 * contains all of the functions need to compute and communicate cell gradients.
 * aggregates the functors that are needed to compute and sum boundary and internal contributions
 * also communicates the cell gradient data on ghosted cells.
 */
class GreenGauss {
  typedef typename ViewTypes::scalar_field_type scalar_field_type;
  typedef typename ViewTypes::solution_field_type solution_field_type;
  typedef typename ViewTypes::gradient_field_type gradient_field_type;
  public:
    GreenGauss(){}

    GreenGauss(Faces * internal_faces, std::vector<Faces *> * bc_faces, Cells * cells, struct MeshData * mesh_data, int total_send_count, int total_recv_count):
      internal_faces_(internal_faces),
      bc_faces_(bc_faces),
      cells_(cells),
      mesh_data_(mesh_data)
      // ghosted_gradient_vars("ghosted_gradient_vars", total_recv_count*5*3),
      // ghosted_gradient_vars_host(ghosted_gradient_vars),
      // shared_gradient_vars("shared_gradient_vars", total_send_count*5*3),
      // shared_gradient_vars_host(shared_gradient_vars)
        {}

    //computes the gradient on locally owned cells.
    void compute_gradients(solution_field_type sol_np1_vec, gradient_field_type gradients){
      printf("Gradients (1)\n");

      //Internal Faces
      const int ninternal_faces = internal_faces_->nfaces_;
      green_gauss_face face_gradient(*internal_faces_, sol_np1_vec, *cells_);
      for(int i = 0; i < ninternal_faces; i++) face_gradient(i);
      // Kokkos::parallel_for(ninternal_faces, face_gradient);

      printf("Gradients (2)\n");

      //Boundary Faces
      typename std::vector<Faces *>::iterator bcf_iter, bcf_iter_end;
      bcf_iter = bc_faces_->begin();
      bcf_iter_end = bc_faces_->end();
      for(; bcf_iter != bcf_iter_end; ++bcf_iter){
        Faces * faces = *bcf_iter;
        const int nboundary_faces = faces->nfaces_;
        green_gauss_boundary_face bc_gradient(*faces, sol_np1_vec, *cells_);
        for(int i = 0; i < nboundary_faces; i++) bc_gradient(i);
        // Kokkos::parallel_for(nboundary_faces, bc_gradient);
      }

      printf("Gradients (3)\n");

      //Sum of all contributions.
      green_gauss_gradient_sum gradient_sum(*cells_, gradients);
      for(int i = 0; i < mesh_data_->num_owned_cells; i++)
        gradient_sum(i);

      printf("Gradients (4)\n");
      // Kokkos::parallel_for(mesh_data_->num_owned_cells, gradient_sum);
      // Kokkos::fence();
    }

    //communicate the computed gradient for ghost cells.
    /*void communicate_gradients(gradient_field_type gradients){
      //copy values to be send from device to host
      extract_shared_tensor<5, 3> extract_shared_gradients(gradients, mesh_data_->send_local_ids, shared_gradient_vars);//sol_np1_vec, send_local_ids, shared_cells);
      // Kokkos::parallel_for(mesh_data_->num_ghosts,extract_shared_gradients);
      // Kokkos::fence();
      // Kokkos::deep_copy(shared_gradient_vars_host, shared_gradient_vars);

      communicate_ghosted_cell_data(mesh_data_->sendCount, mesh_data_->recvCount, shared_gradient_vars_host.ptr_on_device(),ghosted_gradient_vars_host.ptr_on_device(), 15);

      //copy values to be sent from host to device
      // Kokkos::deep_copy(ghosted_gradient_vars, ghosted_gradient_vars_host);
      // insert_ghost_tensor<5, 3> insert_ghost_gradients(gradients, mesh_data_->recv_local_ids, ghosted_gradient_vars);
      // Kokkos::parallel_for(mesh_data_->num_ghosts, insert_ghost_gradients);
      // Kokkos::fence();
    }*/

  private:
    Faces * internal_faces_;
    std::vector<Faces *> * bc_faces_;
    Cells * cells_;
    struct MeshData * mesh_data_;
    // scalar_field_type ghosted_gradient_vars;
    // typename scalar_field_type::HostMirror ghosted_gradient_vars_host;
    // scalar_field_type shared_gradient_vars;
    // typename scalar_field_type::HostMirror shared_gradient_vars_host;
};

#endif