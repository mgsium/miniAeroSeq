#ifndef INCLUDE_MAIN_H_
#define INCLUDE_MAIN_H_

#include <vector>
#include "Options.h"

// Class forward declarations
struct MeshData;

//Function to run on host device.
void run(const Options & simulation_options);

#endif