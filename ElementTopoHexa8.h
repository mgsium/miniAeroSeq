#ifndef INCLUDE_ELEMENT_TOPO_HEXA8_H_
#define INCLUDE_ELEMENT_TOPO_HEXA8_H_

#include "ElementTopo.h"

/// 2x2x2 integration for hexa elements
static const double hex_xi_2x2x2[] =
{ -0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626, -0.577350269189626,
      0.577350269189626, 0.577350269189626, -0.577350269189626 };
static const double hex_eta_2x2x2[] =
{ -0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626,
      -0.577350269189626, 0.577350269189626, 0.577350269189626 };
static const double hex_zeta_2x2x2[] =
{ -0.577350269189626, -0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626,
      0.577350269189626, 0.577350269189626, 0.577350269189626 };
static const double hex_wts_2x2x2[] =
{ 1.0000000000000000, 1.000000000000000, 1.0000000000000000, 1.0000000000000000, 1.0000000000000000,
      1.000000000000000, 1.0000000000000000, 1.0000000000000000 };


/* ElemenTopoHexa
 * class containing the element information for a Hex8 element.
 */
class ElementTopoHexa8: public ElementTopo
{
private:

protected:
   static const double* xi_;
   static const double* eta_;
   static const double* zeta_;

public:
   ElementTopoHexa8();

   virtual ~ElementTopoHexa8();

   unsigned GetNumNodes() const
   {
      return 8;
   }


   int Eval_detJ(const unsigned& gauss_pt, const double* ex, const double* ey, const double* ez, double& detJ,
         double* J, double* dNdxi, double* dNdeta, double* dNdzeta, double* dummy) const;

   void Eval_dNdn(const unsigned& gauss_pt, double* dNdxi, double* dNdeta, double* dNdzeta,
         double* dummy) const;

   double ComputeVolume(const double* ex, const double* ey, const double* ez) const;
};

#endif
