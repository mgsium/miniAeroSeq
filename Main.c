#include "Main.h"
#include "TimeSolverExplicitRK4.h"
#include "Options.h"
#include "MeshData.h"

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
    // Parallel3DMesh mesh(nx, ny, nz, lx, ly, lz, problem_type, angle);
    struct MeshData mesh_data;
    // mesh.fillMeshData<Kokkos::DefaultExecutionSpace>(mesh_data);
    TimeSolverExplicitRK4 * time_solver = new TimeSolverExplicitRK4(mesh_data, simulation_options);
    time_solver->Solve();
    delete time_solver;
}