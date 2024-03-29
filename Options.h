#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <fstream>
#include <iostream>

/*Options
 * struct containing options that are read-in from the input file.
 * The input file is always named miniaero.inp
 */

struct Options{

  int problem_type;
  double lx, ly, lz, angle;
  int nx, ny, nz;
  int ntimesteps;
  double dt;
  int output_results;
  int output_frequency;
  int second_order_space;
  int viscous;

  Options():
    problem_type(0),
    nx(10),
    ny(10),
    nz(10),
    ntimesteps(1),
    dt(5e-8),
    output_results(0),
    output_frequency(10),
    second_order_space(0)
  {

  }

  void read_options_file(){
    //Options file should look like:
    // problem_type (0 - Sod, 1 - Viscous Flat Plate, 2 - Inviscid Ramp)
    // lx ly lz ramp_angle (either SOD(angle=0)  or ramp problem)
    // nx ny nz
    // ntimesteps
    // dt
    // nthreads
    // output_results (0 - no, anything else yes)
    // information output_frequency
    // second order space (0 - no, anything else yes)
    // viscous (0 - no, anything else yes)

    std::ifstream option_file( "miniaero.inp" );

    if(option_file){}
    else {std::cout << "miniaero.inp does not exist." << std::endl;}

    option_file >> problem_type;
    option_file >> lx >> ly >> lz >> angle;
    option_file >> nx >> ny >> nz;
    option_file >> ntimesteps;
    option_file >> dt;
    option_file >> output_results;
    option_file >> output_frequency;
    option_file >> second_order_space;
    option_file >> viscous;
    option_file.close();
  }
};

#endif