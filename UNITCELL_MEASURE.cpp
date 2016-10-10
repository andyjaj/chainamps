/**
 *@file UNITCELL_MEASURE.cpp Makes measurements on a translationally invariant system.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "ajaj_common.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "measurement.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){


  ajaj::QNVector charge_rules;
  ajaj::Basis basis;
  ajaj::UnitCell old(basis);

  ajaj::iMEAS_Args RuntimeArgs(argc,argv);
  for (auto&& f : RuntimeArgs.files()){


    std::ifstream infile;
    infile.open(f.c_str(),ios::in | ios::binary);
    if (infile.is_open()){
      std::cout << "Loading " << f << std::endl;

      ajaj::UnitCell AA(ajaj::load_UnitCell_binary(infile,charge_rules,basis));

      for (auto&& m : AA.Matrices){
	m.print_indices();
	m.print_sparse_info();

      }
      
      for (auto&& l : AA.Lambdas){
	std::cout << ajaj::entropy(l) <<std::endl;
      }
      if (old.size()!=0){
	std::cout << ajaj::Overlap(old,AA) <<std::endl;
      }
      old=std::move(AA);

    }
    else {
      std::cout << "Couldn't open " << f << std::endl;
      return 0;
    }
    
  }
  return 1;
}
