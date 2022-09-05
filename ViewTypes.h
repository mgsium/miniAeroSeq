#ifndef INCLUDE_VIEWTYPES_H_
#define INCLUDE_VIEWTYPES_H_

#include <memory>

/*ViewTypes
 * struct that contains typedefs of commonly used Kokkos Views.
 */
struct ViewTypes{
  typedef double * scalar_field_type;
  typedef double * solution_field_type[5];
  typedef int * face_cell_conn_type[2];
  typedef int ** cell_face_conn_type;
  typedef double * vector_field_type[3];
  typedef double ** cell_storage_field_type[5];
  typedef double * gradient_field_type[5][3];
  typedef double ** gradient_storage_field_type[5][3];

  typedef const double * c_scalar_field_type;
  typedef const double * c_solution_field_type[5];
  typedef const int * c_face_cell_conn_type[2];
  typedef const int ** c_cell_face_conn_type;
  typedef const double * c_vector_field_type[3];
  typedef const double ** c_cell_storage_field_type[5];
  typedef const double * c_gradient_field_type[5][3];
  typedef const double ** c_gradient_storage_field_type[5][3];

  typedef const double * c_rnd_scalar_field_type;
  typedef const double * c_rnd_solution_field_type[5];
  typedef const int * c_rnd_face_cell_conn_type[2]; 
  typedef const int ** c_rnd_cell_face_conn_type;
  typedef const double *  c_rnd_vector_field_type[3];
  typedef const double ** c_rnd_cell_storage_field_type[5];
  typedef const double *  c_rnd_gradient_field_type[5][3];
  typedef const double ** c_rnd_gradient_storage_field_type[5][3];
};

#endif
