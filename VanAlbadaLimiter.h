#ifndef _INCLUDE_VANALBADALIMITER_H_
#define _INCLUDE_VANALBADALIMITER_H_

#include <float.h>
#include "MathToolsDevice.h"

/* VanAlbadaLimiter
 * Class with single function that computes the the VanAlbada limiter value.
 */
class VanAlbadaLimiter
{

public:

 static double limit(double dumax, double dumin, double du){
 double yval = 2;

  if (du > DBL_EPSILON)
  {
    yval = dumax/du;
  }
  else if (du < -DBL_EPSILON)
  {
    yval = dumin/du;
  }

  double phi = 1;
  if (yval < 2){
    phi = (4*yval-yval*yval)/(yval*yval-4*yval+8);
    phi = MathTools::max(phi,0.0);
    phi = MathTools::min(phi,1.0);
  }

  return phi;
  }
};

#endif
