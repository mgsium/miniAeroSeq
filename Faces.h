#ifndef INCLUDE_FACES_H_
#define INCLUDE_FACES_H_

#include "Face.h"
#include "ViewTypes.h"

/*Faces
 * struct containing the face data used for simulation.
 * Includes things such as coordinates, face_normal, and face to element
 * connectivity.
 */
struct Faces{
   typedef int * id_type;
   typedef typename ViewTypes::vector_field_type vector_field_type;
   typedef typename ViewTypes::face_cell_conn_type face_cell_conn_type;

public:
   int nfaces_;
   int ncells_;
   vector_field_type coordinates_, face_normal_, face_tangent_, face_binormal_;
   face_cell_conn_type face_cell_conn_;
   face_cell_conn_type cell_flux_index_;
   id_type permute_vector_;

   Faces(){}

  // Fix memory allocation on initialization
   Faces(int nfaces, int ncells) :
   nfaces_(nfaces),
   ncells_(ncells)
   /*coordinates_(nfaces),
   face_normal_(nfaces),
   face_tangent_(nfaces),
   face_binormal_(nfaces),
   face_cell_conn_(nfaces),
   cell_flux_index_(nfaces),
   permute_vector_(nfaces)*/
  {
    for(int i=0; i < 3; i++) {
      coordinates_[i] = new double[nfaces];
      face_normal_[i] = new double[nfaces];
      face_tangent_[i] = new double[nfaces];
      face_binormal_[i] = new double[nfaces];
    }
    for(int i=0; i < 2; i++) {
      face_cell_conn_[i] = new int[nfaces];
      cell_flux_index_[i] = new int[nfaces];
    }
    permute_vector_ = new int[nfaces];
  }

};

/* copy_faces
 * functor to copy from host setup datastructure
 * to Kokkos datastructures.
 */
/*void copy_faces(Faces device_faces, std::vector<Face> & mesh_faces){

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
  }
}*/

#endif