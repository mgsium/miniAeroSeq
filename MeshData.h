#ifndef INCLUDE_MESH_DATA_H_
#define INCLUDE_MESH_DATA_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Cells.h"
#include "Faces.h"
// struct Cells;
// struct Faces;

/*MeshData
 * Struct that contains needed mesh data
 * Includes ghosting information.
 */
struct MeshData {
  int num_ghosts;
  int num_owned_cells;

  // Host data
  std::vector<int> send_offset, recv_offset;
  std::vector<int> sendCount, recvCount;

  // Device data
  typedef int* id_map_type;
  id_map_type send_local_ids, recv_local_ids;
  Cells mesh_cells;
  Faces internal_faces;
  std::vector<std::pair<std::string, Faces> > boundary_faces;
};

#endif
