#include "Face.h"
#include "MathTools.h"

#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>

Face::Face(const int cv_minus, const int cv_plus, int elem1_flux_index,
        int elem2_flux_index, const int& num_face_nodes,
        std::vector<int> face_nodes, std::vector<double> & node_coordinates) :
      elem1_(cv_minus), elem2_(cv_plus), elem1_flux_index_(elem1_flux_index), elem2_flux_index_(elem2_flux_index), area_(0.0), face_nodes_(face_nodes)
{
   coords_[0] = 0.0;
   coords_[1] = 0.0;
   coords_[2] = 0.0;

   const int len = num_face_nodes;
   double n_xyz[len][3];

   for (int i = 0; i < num_face_nodes; ++i)
   {
       for(int j = 0; j < 3; ++j){
       n_xyz[i][j] = node_coordinates[face_nodes_[i]*3+j];
       }
   }


   for (int in = 0; in < num_face_nodes; ++in)
   {
       MathTools::Vec3SumInto(n_xyz[in], coords_);
   }
   MathTools::Vec3Scale(1.0 / num_face_nodes, coords_);

   const double vec1[] =
   { n_xyz[1][0] - n_xyz[0][0], n_xyz[1][1] - n_xyz[0][1], n_xyz[1][2] - n_xyz[0][2] };
   const double vec2[] =
   { n_xyz[2][0] - n_xyz[0][0], n_xyz[2][1] - n_xyz[0][1], n_xyz[2][2] - n_xyz[0][2] };
   const double vec3[] =
   { n_xyz[3][0] - n_xyz[0][0], n_xyz[3][1] - n_xyz[0][1], n_xyz[3][2] - n_xyz[0][2] };


   double normal1[] =
   { 0.0, 0.0, 0.0 };
   double normal2[] =
   { 0.0, 0.0, 0.0 };

   MathTools::Vec3Cross(vec1, vec2, normal1);
   MathTools::Vec3Cross(vec2, vec3, normal2);
   MathTools::Vec3Average(normal1, normal2, a_vec_);

   //Compute the tangent vector
   double abs_a_vec[] =
   { std::abs(a_vec_[0]), std::abs(a_vec_[1]), std::abs(a_vec_[2]) };
   int i1 = std::max_element(abs_a_vec, abs_a_vec + 3) - abs_a_vec;
   int i2 = i1 + 1;
   int i3 = i1 + 2;
   i2 = (i2 > 2) ? i2 - 3 : i2;
   i3 = (i3 > 2) ? i3 - 3 : i3;

   const double denom = std::sqrt(a_vec_[i1] * a_vec_[i1] + a_vec_[i3] * a_vec_[i3]);

   t_vec_[i2] = 0.0;
   t_vec_[i1] = a_vec_[i3] / denom;
   t_vec_[i3] = -a_vec_[i1] / denom;

   //Compute binormal
   MathTools::Vec3Cross(a_vec_, t_vec_, b_vec_);

}

// ------------------------------------------------------------------------------------------------

Face::~Face()
{
}