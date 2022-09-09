#include "Copy.h"

/*copy_cell_data
 * functor to copy the Cell information from the setup datastructure to
 * Kokkos datastructure.
 */
void copy_cell_data(Cells device_cells, std::vector<Cell> & mesh_cells){
  
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
}

/* copy_faces
 * functor to copy from host setup datastructure
 * to Kokkos datastructures.
 */
void copy_faces(Faces device_faces, std::vector<Face> & mesh_faces){

  //Need host mirror
  typedef int * id_type;
  typedef typename ViewTypes::vector_field_type vector_field_type;
  typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;

  int nfaces = mesh_faces.size();

    double a_vec[3], t_vec[3], b_vec[3];
    for(int i = 0; i < mesh_faces.size(); ++i){
        device_faces.face_cell_conn_[i][0]  = mesh_faces[i].GetElem1();
        device_faces.face_cell_conn_[i][1]  = mesh_faces[i].GetElem2();
        device_faces.cell_flux_index_[i][0] = mesh_faces[i].GetElem1_FluxIndex();
        device_faces.cell_flux_index_[i][1] = mesh_faces[i].GetElem2_FluxIndex();

        const double * coords = mesh_faces[i].GetCoords();
        mesh_faces[i].GetAreaVecAndTangentVec(a_vec, t_vec, b_vec);
        for(int j = 0; j < 3; ++j){
            device_faces.coordinates_[i][j]   = coords[j];
            device_faces.face_normal_[i][j]   = a_vec[j];
            device_faces.face_tangent_[i][j]  = t_vec[j];
            device_faces.face_binormal_[i][j] = b_vec[j];
        }
    }

  /*if(device_faces.face_cell_conn_.extent(0) > 0) {
    typedef Kokkos::View<int *, Kokkos::LayoutStride, Device> view_type;
    typedef Kokkos::BinOp1D< view_type > CompType;
    view_type face_cell_left = Kokkos::subview(device_faces.face_cell_conn_,Kokkos::ALL(),0);

    typedef Kokkos::MinMax<int,Device> reducer_type;
    typedef typename reducer_type::value_type minmax_type;
    minmax_type minmax;
    Kokkos::parallel_reduce(face_cell_left.extent(0), KOKKOS_LAMBDA (const int& i, minmax_type& lminmax) {
      if(face_cell_left(i)<lminmax.min_val) lminmax.min_val = face_cell_left(i);
      if(face_cell_left(i)>lminmax.max_val) lminmax.max_val = face_cell_left(i);
    },reducer_type(minmax));

    Kokkos::BinSort<view_type, CompType, Device, int> bin_sort(face_cell_left,CompType(face_cell_left.extent(0)/2,minmax.min_val,minmax.max_val),true);
    bin_sort.create_permute_vector();
    Kokkos::deep_copy(device_faces.permute_vector_, bin_sort.sort_order);
  }*/
}