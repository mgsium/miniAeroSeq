#ifndef INCLUDE_FACE_H_
#define INCLUDE_FACE_H_

#include <vector>


class Face
{
private:

protected:

   int elem1_;
   int elem2_;
   int elem1_flux_index_;
   int elem2_flux_index_;
   double coords_[3];
   double a_vec_[3];
   double t_vec_[3];
   double b_vec_[3];


   double area_;

   std::vector<int> face_nodes_;

public:
   Face(const int cv_minus, const int cv_plus, int elem1_flux_index, int elem2_flux_index, const int& num_face_nodes,
           std::vector<int> face_nodes, std::vector<double> & node_coordinates);

   virtual ~Face();

   int GetElem1() const
   {
      return elem1_;
   }
   int GetElem2() const
   {
      return elem2_;
   }

   int GetElem1_FluxIndex() const
   {
      return elem1_flux_index_;
   }
   int GetElem2_FluxIndex() const
   {
      return elem2_flux_index_;
   }

   const std::vector<int> & GetFaceNodes() const
   {
      return face_nodes_;
   }

   const double* GetCoords() const
   {
      return coords_;
   }

   const double* GetAreaVec() const
   {
      return a_vec_;
   }

   void GetAreaVecAndTangentVec(double a_vec[], double t_vec[], double b_vec[]) const
   {
      a_vec[0] = a_vec_[0];
      a_vec[1] = a_vec_[1];
      a_vec[2] = a_vec_[2];

      t_vec[0] = t_vec_[0];
      t_vec[1] = t_vec_[1];
      t_vec[2] = t_vec_[2];

      b_vec[0] = b_vec_[0];
      b_vec[1] = b_vec_[1];
      b_vec[2] = b_vec_[2];
   }
};

#endif
