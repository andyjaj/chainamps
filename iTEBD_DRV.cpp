/**
 *@file iTEBD_DRV.cpp Driver file for iTEBD.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>

#include "ajaj_common.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "TEBD_routines.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){

  ajaj::iTEBD_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    myModel.basis().print();
    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double minS(1.0e-14);
    ajaj::uMPXInt number_of_time_steps(RuntimeArgs.number_of_steps());
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order());
    ajaj::uMPXInt measurement_interval(RuntimeArgs.measurement_interval());
    double time_step(RuntimeArgs.step_size());
    std::vector<ajaj::MPO_matrix> measured_operators;// can insert runtime measurements here

    ajaj::UnitCell Initial(myModel.basis());
    std::string Name;
    std::stringstream Rss;
    Rss<<RuntimeArgs.filename();
    if (!RuntimeArgs.initial_unit_cell().empty()){
      std::ifstream infile;
      infile.open(RuntimeArgs.initial_unit_cell().c_str(),ios::in | ios::binary);
      if (infile.is_open()){
	std::cout << "Using Initial UnitCell " << RuntimeArgs.initial_unit_cell() <<std::endl;
	Initial=ajaj::load_UnitCell_binary(infile,myModel.basis().getChargeRules(),myModel.basis());//populates basis
	std::istringstream fnss(RuntimeArgs.initial_unit_cell());
	std::string filename;
	getline(fnss,filename,'.');
	Name=ajaj::StripName(filename);
	Rss<<"_"<<Name;
      }
      else {
	std::cout << "Error: couldn't open file " << RuntimeArgs.initial_unit_cell() <<std::endl;
	return 1;
      }
    }
    else {
      Name="Ortho";
      Initial=MakeProductStateUnitCell(myModel.basis(),0,ajaj::State(myModel.basis().getChargeRules()));
      //Initial=MakeProductStateUnitCell(myModel.basis(),std::vector<std::pair<ajaj::uMPXInt,double> > ({{0,1.0/sqrt(2.0)},{1,1.0/sqrt(2.0)}}),ajaj::State(myModel.basis().getChargeRules()));
    }


    ajaj::DataOutput results(ajaj::OutputName(Rss.str(),"Evolution.dat"),"Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap)");

      //ajaj::UnitCell Initial(MakeProductStateUnitCell(myModel.basis(),0,ajaj::State(myModel.basis().getChargeRules())));
    //const ajaj::UnitCell Initial(MakeProductStateUnitCell(myModel.basis(),std::vector<std::pair<ajaj::uMPXInt,double> > ({{0,1.0/sqrt(2.0)},{1,1.0/sqrt(2.0)}}),ajaj::State(myModel.basis().getChargeRules())));

    ajaj::iTEBD infrun(myModel.H_MPO,Initial,time_step,results,Name,trotter_order);
    infrun.evolve(number_of_time_steps,measured_operators,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
    return 0;
  }
  else {
    return 1;
  }
}

