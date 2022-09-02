#ifndef INCLUDE_TIMESOLVER_EXPLICIT_RK4_H_
#define INCLUDE_TIMESOLVER_EXPLICIT_RK4_H_

// C++ system files
#include <cstdio>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>

// TPL header files
#include "MemoryUsage.h"

//Input File Options
#include "Options.h"

//Data structure headers
#include "MeshData.h"
#include "Faces.h"
#include "Cells.h"

// Functor headers
#include "Extrapolate_BC.h"
#include "Tangent_BC.h"
#include "Inflow_BC.h"
#include "NoSlip_BC.h"
#include "Initial_Conditions.h"
#include "Flux.h"
#include "Roe_Flux.h"
#include "Viscous_Flux.h"
// #include "CopyGhost.h"
#include "GreenGauss.h"
#include "StencilLimiter.h"

/* TimeSolverData
 * Class that contains options for time marching.
 */
struct TimeSolverData
{
   unsigned time_it;
   unsigned max_its;
   double start_time;
   double time;
   double dt;

   TimeSolverData() :
      time_it(0),
      max_its(1),
      start_time(0.0),
      time(0.0),
      dt(5e-8)
   {

   }

   ~TimeSolverData()
   {

   }
};

/* update
 * functor that updates solution using old solution, residual and scaling of residual.
 */
void update(
  double alpha, 
  ViewTypes::solution_field_type res, 
  ViewTypes::solution_field_type soln,
  ViewTypes::solution_field_type solnp1, 
  int i
) {
    for(int icomp=0; icomp<5; icomp++)
    {
      solnp1[i][icomp] = soln[i][icomp] + alpha * res[i][icomp];
    }
}

/* copy
 * functor that copies from one solution array to another.
 */
void copy (
  ViewTypes::solution_field_type soln_src, 
  ViewTypes::solution_field_type soln_dst,
  int i
) {
  for(int icomp=0; icomp<5; icomp++)
  {
    soln_dst[i][icomp] = soln_src[i][icomp];
  }
}

/*TimeSolverExplicitRK4
 * Class that runs RK4 solver which is the main loop in this code
 * All of the physics kernels are called from the inner RK4 stage loop.
 * Is a basic implementation of the standard RK4 solver.
 */
class TimeSolverExplicitRK4
{
public:
   TimeSolverExplicitRK4(struct MeshData & input_mesh_data, const Options & options);
   ~TimeSolverExplicitRK4();
   void Solve();
   bool CheckStopCriteria();

private:
  struct MeshData & mesh_data_;
  TimeSolverData ts_data_;
  unsigned stages_;
  double alpha_[4];
  double beta_[4];
  Options options_;
  TimeSolverExplicitRK4();
};

TimeSolverExplicitRK4::TimeSolverExplicitRK4(struct MeshData & input_mesh_data, const Options & options):
  ts_data_(),
  mesh_data_(input_mesh_data),
  options_(options)
{
  ts_data_.max_its=options_.ntimesteps;
  ts_data_.dt=options_.dt;
  stages_=4;
  alpha_[0] = 0.0;
  alpha_[1] = 1.0/2.0;
  alpha_[2] = 1.0/2.0;
  alpha_[3] = 1.0;
  beta_[0] = 1.0/6.0;
  beta_[1] = 1.0/3.0;
  beta_[2] = 1.0/3.0;
  beta_[3] = 1.0/6.0;
}

TimeSolverExplicitRK4::~TimeSolverExplicitRK4()
{

}

///////////////////////////////////////////////////////////////////////////////////////////////////
void TimeSolverExplicitRK4::Solve()
{
    typedef typename ViewTypes::scalar_field_type scalar_field_type;
    typedef typename ViewTypes::solution_field_type solution_field_type;
    typedef typename ViewTypes::vector_field_type vector_field_type;
    typedef typename ViewTypes::gradient_field_type gradient_field_type;
    typedef int * id_map_type;

    const double midx=options_.lx/2.0;
    double inflow_state[5];
    inflow_state[0]=0.5805;
    inflow_state[1]=503.96;
    inflow_state[2]=0.0;
    inflow_state[3]=0.0;
    inflow_state[4]=343750.0;

    //Faces - Interior and BC
    //Internal Faces
    Faces internal_faces = mesh_data_.internal_faces;
    const int ninternal_faces= internal_faces.nfaces_;

    //Boundary Faces
    std::vector<Faces > extrapolate_faces, tangent_faces, inflow_faces, noslip_faces;
    std::vector<Faces * > bc_faces;

    typename std::vector<std::pair<std::string, Faces > >::iterator bc_iter, bc_iter_end;
    bc_iter = mesh_data_.boundary_faces.begin();
    bc_iter_end = mesh_data_.boundary_faces.end();
    for(;bc_iter != bc_iter_end; bc_iter++){
        bc_faces.push_back(&(bc_iter->second));
        if (bc_iter->first == "Extrapolate")
            extrapolate_faces.push_back(bc_iter->second);
        else if(bc_iter->first == "Tangent")
            tangent_faces.push_back(bc_iter->second);
        else if(bc_iter->first == "Inflow")
            inflow_faces.push_back(bc_iter->second);
        else if(bc_iter->first == "NoSlip")
            noslip_faces.push_back(bc_iter->second);
    }

    //Cells
    Cells cells = mesh_data_.mesh_cells;
    const int nowned_cells = mesh_data_.num_owned_cells;
    const int num_ghosts = mesh_data_.num_ghosts;
    const int ncells = cells.ncells_;

    //Solution Variables
    solution_field_type res_vec("residual", ncells);
    ViewTypes::solution_field_type sol_n_vec("solution_n", ncells);
    solution_field_type sol_np1_vec("solution_np1", ncells);
    solution_field_type sol_temp_vec("solution_temp", ncells); //Needed for RK4 Stages.

    gradient_field_type gradients;
    solution_field_type limiters;
    //Careful not to allocate the views unless they are needed.
    if(options_.second_order_space || options_.viscous){
        gradients = gradient_field_type("gradients", ncells);
    }
    if(options_.second_order_space){
        limiters = solution_field_type("limiters", ncells);
    }

    solution_field_type solution_vec   = sol_n_vec;
    solution_field_type residuals_host = res_vec;
    vector_field_type coordinates_host = cells.coordinates_;

    //setup ghosting information
    int num_procs_, my_id_;
    int total_send_count=0, total_recv_count=0;

    // Gradients
    GreenGauss green_gauss_gradient;

    // Stencil Limiter
    StencilLimiter stencil_limiter;

   if(options_.second_order_space || options_.viscous)
   {
      green_gauss_gradient = GreenGauss(&internal_faces, &bc_faces, &cells, &mesh_data_, total_send_count, total_recv_count);
   }
   if(options_.second_order_space)
   {
      stencil_limiter = StencilLimiter(&internal_faces, &bc_faces, &cells, &mesh_data_, total_send_count, total_recv_count);
   }

   if(options_.problem_type == 0)
   {
      initialize_sod3d(cells, sol_n_vec, sol_temp_vec, midx);
      // Kokkos::parallel_for(nowned_cells, init_fields);
    }
   else
    {
      initialize_constant(cells, sol_n_vec, sol_temp_vec, &inflow_state[0]);
      // Kokkos::parallel_for(nowned_cells, init_fields);
    }
    // Initialize the value for np1 solution which will be updated each RK stage and then copied to n solution at end of timestep.
    // copy<Device> copy_solution( sol_n_vec, sol_np1_vec);
    // Kokkos::parallel_for(nowned_cells, copy_solution);
    // Kokkos::fence();

    for (ts_data_.time_it = 1; ts_data_.time_it <= ts_data_.max_its; ++ts_data_.time_it)
   {
      // Increment the time, do not need to worry about updating it for stages because no source terms depend on the time.
      ts_data_.time += ts_data_.dt;

      // Print time step info
      /*if(ts_data_.time_it % options_.output_frequency == 0  && my_id_==0){
        fprintf(stdout, "\nTime Step #%i:  Time = %16.9e; dt = %16.9e\n", ts_data_.time_it, ts_data_.time,
             ts_data_.dt);
      }*/

      // R-K stages loop ****************************************************
      for (unsigned irk = 0; irk < stages_; ++irk)
      {
        //Update temporary solution used to evaluate the residual for this RK stage
        for (int i=0; i < nowned_cells; i++)
          update(alpha_[irk], res_vec, sol_n_vec, sol_temp_vec, i);

        //Zero fluxes
        zero_cell_flux zero_flux(cells);
        // Kokkos::parallel_for(nowned_cells, zero_flux);
        // Kokkos::fence();

        //Compute Gradients and Limiters
        if(options_.second_order_space || options_.viscous){
          green_gauss_gradient.compute_gradients(sol_temp_vec, gradients);
        }

        if(options_.second_order_space){
          stencil_limiter.compute_min_max(sol_temp_vec);
        }
      }

      //Compute internal face fluxes
      roe_flux inviscid_flux_evaluator;
      if(options_.viscous){
        newtonian_viscous_flux viscous_flux_evaluator;
        if(options_.second_order_space){
          compute_face_flux<true, roe_flux, newtonian_viscous_flux >(internal_faces, sol_temp_vec, gradients, limiters, cells, inviscid_flux_evaluator, viscous_flux_evaluator);
          // Kokkos::parallel_for(ninternal_faces,fluxop);
        }
        else{
          compute_face_flux<false, roe_flux, newtonian_viscous_flux>(internal_faces, sol_temp_vec, gradients, limiters, cells, inviscid_flux_evaluator, viscous_flux_evaluator);
          // Kokkos::parallel_for(ninternal_faces,fluxop);
        }
        // Kokkos::fence();
      } 
      else{
        no_viscous_flux viscous_flux_evaluator;
        if(options_.second_order_space){
          compute_face_flux<true, roe_flux, no_viscous_flux >(internal_faces, sol_temp_vec, gradients, limiters, cells, inviscid_flux_evaluator, viscous_flux_evaluator);
          // Kokkos::parallel_for(ninternal_faces,fluxop);
        }
        else{
          compute_face_flux<false, roe_flux, no_viscous_flux >(internal_faces, sol_temp_vec, gradients, limiters, cells, inviscid_flux_evaluator, viscous_flux_evaluator);
          // Kokkos::parallel_for(ninternal_faces,fluxop);
        }
        // Kokkos::fence();
      }

      //Extrapolated BC fluxes
      typename std::vector<Faces>::iterator ef_iter, ef_iter_end;
      ef_iter = extrapolate_faces.begin();
      ef_iter_end = extrapolate_faces.end();
      for(; ef_iter != ef_iter_end; ++ef_iter){
        Faces bc_faces = *ef_iter;
        const int nboundary_faces = bc_faces.nfaces_;
        compute_extrapolateBC_flux<roe_flux >(bc_faces, sol_temp_vec, cells, inviscid_flux_evaluator);
        // Kokkos::parallel_for(nboundary_faces,boundary_fluxop);
      }
      // Kokkos::fence();

      //Tangent BC fluxes
      typename std::vector<Faces>::iterator tf_iter, tf_iter_end;
      tf_iter = tangent_faces.begin();
      tf_iter_end = tangent_faces.end();
      for(; tf_iter != tf_iter_end; ++tf_iter){
        Faces bc_faces = *tf_iter;
        const int nboundary_faces = bc_faces.nfaces_;
        compute_tangentBC_flux<roe_flux >(bc_faces, sol_temp_vec, cells, inviscid_flux_evaluator);
        // Kokkos::parallel_for(nboundary_faces,boundary_fluxop);
      }
      // Kokkos::fence();

      //Noslip BC fluxes
      typename std::vector<Faces>::iterator if_iter, if_iter_end;
      if_iter = noslip_faces.begin();
      if_iter_end = noslip_faces.end();
      for(; if_iter != if_iter_end; ++if_iter){
        newtonian_viscous_flux viscous_flux_evaluator;
        Faces bc_faces = *if_iter;
        const int nboundary_faces = bc_faces.nfaces_;
        compute_NoSlipBC_flux<roe_flux, newtonian_viscous_flux > boundary_fluxop(bc_faces, sol_temp_vec, cells, inviscid_flux_evaluator, viscous_flux_evaluator);
        // Kokkos::parallel_for(nboundary_faces,boundary_fluxop);
      }
      // Kokkos::fence();

      //Inflow BC fluxes
      typename std::vector<Faces >::iterator nsf_iter, nsf_iter_end;
      nsf_iter = inflow_faces.begin();
      nsf_iter_end = inflow_faces.end();
      for(; nsf_iter != nsf_iter_end; ++nsf_iter){
        Faces bc_faces = *nsf_iter;
        const int nboundary_faces = bc_faces.nfaces_;
        compute_inflowBC_flux<roe_flux > boundary_fluxop(bc_faces, sol_temp_vec, cells, &inflow_state[0], inviscid_flux_evaluator);
        // Kokkos::parallel_for(nboundary_faces,boundary_fluxop);
      }
      // Kokkos::fence();

      //Sum up all of the contributions
      apply_cell_flux flux_residual(cells, res_vec, ts_data_.dt);
      // Kokkos::parallel_for(nowned_cells, flux_residual);
      // Kokkos::fence();

      //Update np1 solution with each stages contribution
      update_fields(beta_[irk],res_vec,sol_np1_vec,sol_np1_vec);
      // Kokkos::parallel_for(nowned_cells, update_fields);
      // Kokkos::fence();
    }

    // Update the solution vector after having run all of the RK stages.
    for (int i=0; i < nowned_cells; i++)
      copy(sol_np1_vec, sol_n_vec, i);
    // Kokkos::parallel_for(nowned_cells, copy_solution);

  /*if(my_id_==0){
     fprintf(stdout,"\n ... Device Run time: %8.2f seconds ...\n", timer.seconds());
   }*/

  size_t current_mem_usage, high_water_mem_usage;
  get_memory_usage(current_mem_usage, high_water_mem_usage);
  current_mem_usage = current_mem_usage/(1024.0*1024.0);
  high_water_mem_usage = high_water_mem_usage/(1024.0*1024.0);
  fprintf(stdout,"\n CPU Memory Usage (Current, High Water) - end of calculation: %lu MB, %lu MB", current_mem_usage, high_water_mem_usage);

  // Output to file on the host.  Requires a deep copy from device to host.
   if(options_.output_results){
     std::ofstream output_file;
     std::stringstream fs;
     fs << "results.";
     // fs << my_id_;
     std::string filename = fs.str();
     output_file.open(filename.c_str(), std::ios::out);

     for(int i=0; i<nowned_cells; i++)
     {
       output_file << cells.coordinates_[i][0] << "\t";
       output_file << cells.coordinates_[i][1] << "\t";
       output_file << cells.coordinates_[i][2] << "\t";
       for(int icomp=0; icomp<5; icomp++){
         output_file << sol_n_vec[i][icomp] << "\t";
       }
       output_file << std::endl;
     }
     output_file.close();
   }
}
