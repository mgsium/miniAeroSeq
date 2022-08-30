#ifndef _INCLUDE_VENKATLIMITER_H_
#define _INCLUDE_VENKATLIMITER_H_

#include <float.h>
#include "MathToolsDevice.h"

/* VenkatLimiter
 * Class with single function that computes the the Venkat limiter value.
 */
class VenkatLimiter
{

public:

 static double limit(double dumax, double dumin, double du, double deltax3){
  double beta = 1;
  double epstilde2 = deltax3*beta*beta*beta;
  double phi = 0;

  double denom = 0;
  double num = 0;


  if (du > 1e-40)
  {
    num = (dumax*dumax + epstilde2)*du + 2*du*du*dumax;
    denom = du*(dumax*dumax  + 2*du*du + dumax*du + epstilde2);
    phi = num/denom;
  }
  else if (du < -1e-40)
  {
    num  = (dumin*dumin + epstilde2)*du + 2*du*du*dumin;
    denom = du*(dumin*dumin  + 2*du*du + dumin*du + epstilde2);
    phi = num/denom;
  }
  else
  {
    phi = 1;
  }


  return phi;
  }
};

#endif
