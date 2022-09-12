#include "Main.h"
#include "TimeSolverExplicitRK4.h"
#include "Parallel3DMesh.h"
#include "MeshData.h"
#include "Options.h"

#include <unistd.h>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

int main(int argc, char *argv[])
{
    Options simulation_options;
    simulation_options.read_options_file();
    run(simulation_options);

    return 0;
}

void run(const Options & simulation_options)
{
    int num_procs, my_id;

    // Setup mesh on host
    int nx = simulation_options.nx, ny = simulation_options.ny, nz = simulation_options.nz;
    double lx = simulation_options.lx, ly = simulation_options.ly, lz = simulation_options.lz;
    double angle = simulation_options.angle;
    int problem_type = simulation_options.problem_type;
    
    // Fill mesh
    Parallel3DMesh mesh(nx, ny, nz, lx, ly, lz, problem_type, angle);
    struct MeshData mesh_data;
    printf("Created Mesh");

    mesh.fillMeshData(mesh_data);
    printf("Owned Cells: %d\n", mesh_data.num_owned_cells);
    printf("Filled Mesh\n");

    TimeSolverExplicitRK4 * time_solver = new TimeSolverExplicitRK4(mesh_data, simulation_options);
    printf("Initialized\n");

    time_solver->Solve();
    printf("Solved\n");

    delete time_solver;
}
