#ifndef INCLUDE_MESH_PROCESSOR_H_
#define INCLUDE_MESH_PROCESSOR_H_
#include <vector>
#include <set>

class Face;
class Cell;

/*FaceData
 * struct containing face data that is used for setup
 * not used for computation.
 */
struct FaceData
{
   int cv_minus;
   int cv_plus;

   int flux_index_minus;
   int flux_index_plus;

   std::vector<int> node_lids;

   int min_node_lid;

   FaceData() :
         cv_minus(-1), cv_plus(-1), flux_index_minus(-1), flux_index_plus(-1)
   {
   }

   ~FaceData()
   {
   }
};


/*Free functions used for mesh setup*/

void create_faces(std::vector<int> & element_node_conn, std::vector<Face> & faces, std::vector<double> & node_coordinates, int num_elems, int num_nodes);

void compute_cell_volumes(std::vector<Cell> & cells, std::vector<int> & element_node_conn, std::vector<double> & node_coordinates, int num_elems, int num_ghosted);

void compute_cell_centroid(std::vector<Cell> & cells, std::vector<int> & element_node_conn, std::vector<double> & node_coordinates, int num_elems);

void delete_ghosted_faces(std::vector<Face> & all_faces, int ghosted_elem_index);

void organize_BC_faces(std::vector<Face> & all_faces, std::vector<Face> & bc_faces, std::set<int> & bc_nodes);

void extract_BC_faces(std::vector<Face> & all_faces, std::vector<Face> & bc_faces);

#endif
