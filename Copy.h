#ifndef COPY_H_DEFINED
#define COPY_H_DEFINED

#include "Faces.h"
#include "Cells.h"

void copy_cell_data(Cells device_cells, std::vector<Cell> & mesh_cells);
void copy_faces(Faces device_faces, std::vector<Face> & mesh_faces);

#endif