#ifndef INCLUDE_CELL_H_
#define INCLUDE_CELL_H_

/*Cell
 * Class to hold volume and coordinate data for mesh cells.
 * Only used temporarily in setting up the mesh.
 */

class Cell
{
private:

protected:
   double coords_[3];
   double volume_;


public:
   Cell();

   virtual ~Cell();

   double* GetCoords()
   {
      return coords_;
   }

   void SetCoords(const double* xyz)
   {
      coords_[0] = xyz[0];
      coords_[1] = xyz[1];
      coords_[2] = xyz[2];
   }

   double& GetVolume()
   {
      return volume_;
   }

};

#endif
