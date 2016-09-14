/**
 *@file ChainAMPS.cpp Driver file for iDMRG.
 */

#include "model.hpp"
#include "make_model.hpp"
#include "command_line_input.hpp"
#include "states.hpp"
#include "data.hpp"
#include "DMRG_routines.hpp"


int main(int argc, char** argv){
  static double convergence_test=-0.0;

  ajaj::iDMRG_Args RuntimeArgs(argc,argv);
  const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
  myModel.basis().print();
  return 0;

}
