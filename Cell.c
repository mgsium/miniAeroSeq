#include "Cell.h"

// ------------------------------------------------------------------------------------------------

Cell::Cell() :
   volume_(0.0)
{
   coords_[0] = 0.0;
   coords_[1] = 0.0;
   coords_[2] = 0.0;
}

// ------------------------------------------------------------------------------------------------

Cell::~Cell()
{
}
