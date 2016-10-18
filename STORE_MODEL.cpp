/**
 *@file STORE_OPERATORS.cpp Form a model and store the vertex operators as (text) based sparse matrices.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "ajaj_common.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "measurement.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){
  ajaj::Store_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){

    const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));

    std::string input_filename(ajaj::StripName(RuntimeArgs.filename()));
    std::ostringstream bss;
    bss << input_filename <<".BASIS";
    std::ofstream basisfile;
    basisfile.open(bss.str().c_str(),ios::out | ios::trunc);
    if (basisfile.is_open()){
      basisfile << "# Index";
      if (myModel.basis().getChargeRules().size()) basisfile << " {";
      for (auto c : myModel.basis().getChargeRules()){
	std::ostringstream css;
	css << " Z";
	if (c>0){
	  css << "_" << c;
	}
	basisfile << css.str(); 
      }
      if (myModel.basis().getChargeRules().size()) basisfile << " }";
      basisfile << std::endl;
      size_t i(0);
      for (auto&& s :myModel.basis()){
	basisfile << i++ << " "  << s <<std::endl;
      }
      basisfile.close();
    }
    else {
      std::cout << "Couldn't open " << bss.str() << " for writing!" << std::endl;
      return 0;
    }
    for (auto&& O : myModel.vertex.Operators){
      std::ostringstream outss;
      outss << input_filename <<"-"<<O.Name <<".SPARSEMATRIX";
      std::string outfilename(outss.str());
      std::replace( outfilename.begin(), outfilename.end(), ' ', '_');
      std::ofstream outfile;
      outfile.open(outfilename.c_str(),ios::out | ios::trunc);
      if (outfile.is_open()){
	O.MatrixElements.fprint(outfile);
	outfile.close();
      }
      else {
	std::cout << "Error opening " << outss.str() << " for writing!" << std::endl;
	return 0;
      }
    }

    return 1;
  }
  return 0;
}
