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

    std::string subname (RuntimeArgs.filename().find("/")!=std::string::npos ?  RuntimeArgs.filename().substr(RuntimeArgs.filename().rfind("/")+1,RuntimeArgs.filename().length()) : RuntimeArgs.filename());

    const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));

    for (auto&& O : myModel.vertex.Operators){



      std::ostringstream outss;
      outss << subname <<"-"<<O.Name <<".SPARSEMATRIX";
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
