#ifndef INCLUDE_MATH_TOOLS_H_
#define INCLUDE_MATH_TOOLS_H_

#include <cmath>

/*MathTools
 * Class that contains some convenience math functions that are
 * templated to run on device.
 */
class MathTools
{

public:

   // ----------------------------------------------------------------------------
   static double Vec3Norm(const double a[])
   {
      return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
   }

   static void MatVec(const int N, const double alpha, const double A[], const double x[], const double beta, double y[])
   {
      for(int i=0; i < N; ++i) {
        y[i] *= beta;
        for(int j=0; j < N; ++j) {
          y[i] += alpha * A[N*i + j] * x[j];
        }
      }
   }

   static void MatVec5(const double alpha, const double A[], const double x[], const double beta, double y[])
   {
      for(int i=0; i < 5; ++i) {
        y[i] *= beta;
        for(int j=0; j < 5; ++j) {
          y[i] += alpha * A[5*i + j] * x[j];
        }
      }
   }

   static double min(double val1, double val2)
   {
     return std::min(val1, val2);
   }
   static double max(double val1, double val2)
   {
     return std::max(val1, val2);
   }
};
#ifdef KOKKOS_HAVE_CUDA
template <>
double MathTools::min(double val1, double val2)
  {
    return fmin(val1, val2);
  } 

template <>
double MathTools::max(double val1, double val2)
  {
    return fmax(val1, val2);
  } 
#endif

#endif
