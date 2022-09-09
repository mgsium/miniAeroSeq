#ifndef INCLUDE_ELEMENT_TOPO_H_
#define INCLUDE_ELEMENT_TOPO_H_

#include <cassert>

static const double DETJ_EPS = 1.0e-300;

/* --------------------------------------------------------------------------------------------- */

/*ElemenTopy
 * base class for element topologies
 * stores and computes things such as number of nodes, gauss points, and volume.
 */

class ElementTopo
{
private:

protected:
   unsigned num_gp_;
   const double* wts_;

public:
   ElementTopo();

   virtual ~ElementTopo();

   virtual unsigned GetNumNodes() const = 0;

   unsigned GetNumGaussPoints() const
   {
      return num_gp_;
   }

   const double* GetGaussWeights() const
   {
      return wts_;
   }


   virtual int Eval_detJ(const unsigned& gauss_pt, const double* ex, const double* ey, const double* ez,
         double& detJ, double* J, double* dNdxi, double* dNdeta, double* dNdzeta, double* dummy) const = 0;

   virtual void Eval_dNdn(const unsigned& gauss_pt, double* dNdn1, double* dNdn2, double* dNdn3,
         double* dNdx4) const = 0;

   virtual double ComputeVolume(const double* ex, const double* ey, const double* ez) const = 0;

};

#endif
