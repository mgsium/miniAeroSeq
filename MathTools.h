#ifndef INCLUDE_MATH_TOOLS_H_
#define INCLUDE_MATH_TOOLS_H_

#include <cmath>
#include <algorithm>



/*MathTools
 * Class that contains some convenience math functions.
 */

class MathTools
{
private:

protected:

public:

   static void Vec3Cross(const double a[], const double b[], double c[])
   {
      c[0] = a[1] * b[2] - b[1] * a[2];
      c[1] = -a[0] * b[2] + b[0] * a[2];
      c[2] = a[0] * b[1] - b[0] * a[1];
   }


   static void Vec3Average(const double a[], const double b[], double c[])
   {
      c[0] = 0.5 * (a[0] + b[0]);
      c[1] = 0.5 * (a[1] + b[1]);
      c[2] = 0.5 * (a[2] + b[2]);
   }


   static void Vec3SumInto(const double a[], double c[])
   {
      c[0] += a[0];
      c[1] += a[1];
      c[2] += a[2];
   }

   static void Vec3Scale(const double scalar, double c[])
   {
      c[0] *= scalar;
      c[1] *= scalar;
      c[2] *= scalar;
   }

   static double VecNDot(const unsigned N, const double a[], const double b[])
   {
      double result = 0;
      for(unsigned i=0; i < N; ++i)
        result += a[i] * b[i];
      return result;
   }


};

#endif
