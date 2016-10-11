/**
 *@file UNITCELL_MEASURE.cpp Makes measurements on a translationally invariant system.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>
#include <algorithm>

#include "ajaj_common.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "sparse_interface.hpp"
#include "measurement.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){

  ajaj::iMEAS_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    //files() returns std::vector<std::string> of files to work on
    //Options:
    //separation() (between operators, indicates measurement at multiple points)
    //operator_filenames() (can be 1 or 2 at the moment, indicates measurement at multiple locations, specified by distance).
    //use_filename_index(), strip and index from filename and use it to order results, default true

    //read in (possible) operators
    std::vector<ajaj::SparseMatrix> MElements;
    if (RuntimeArgs.operator_filenames().size()) {
      for (auto&& opfn : RuntimeArgs.operator_filenames()){
	///open file
	//make into sparse matrix
	//should be coord list, first line dims
	std::ifstream is;
	is.open(opfn.c_str(),std::ios::in);
	if (is.is_open()){
	  std::string fileline;
	  getline(is,fileline);
	  ajaj::MPXInt rows(0);
	  ajaj::MPXInt cols(0);
	  std::istringstream line1ss(fileline);
	  if (line1ss >> rows && line1ss >> cols && line1ss.eof() && rows==cols){
	    MElements.emplace_back(rows,cols);
	    while (getline(is,fileline)){
	      std::istringstream liness(fileline);
	      ajaj::MPXInt row;
	      ajaj::MPXInt col;
	      std::complex<double> value;
	      if (liness >> row && liness >> col && liness >> value){
		if (value!=0.0){
		  MElements.back().entry(row,col,value);
		}
	      }
	      else {
		std::cout << "Malformed operator in " << opfn << std::endl;
		return 0;
	      }
	    }
	    MElements.back().finalise();
	  }
	  else {
	    std::cout << "Malformed dimension info in " << opfn << std::endl;
	    return 0;
	  }
	}
	else {
	  std::cout << "Couldn't open " << opfn << std::endl;
	  return 0;
	}
	//MElements.back().print();
      }
    }

    //check dims
    ajaj::MPXInt dim(MElements.size() ? MElements[0].rows() : 0);

    if (MElements.size()>1){
      for (auto&& m : MElements){
	if (dim!= m.rows()){
	  std::cout << "Operator dimensions don't match! " << dim << " " << m.rows() <<std::endl;
	  return 0;
	}
      }
    }

    ajaj::QNVector charge_rules;
    ajaj::Basis basis;
    ajaj::UnitCell cellbuff(basis);
    std::vector<std::pair<size_t,ajaj::Data> > indexed_results;
    size_t Index(0);
    std::regex ex("_([0-9]+)\\.UNITCELL$"); //regex to match number before filename only.

    for (auto&& f : RuntimeArgs.files()){
      //get an index
      if (RuntimeArgs.use_filename_index()){ //extract from filename
	std::smatch sm;
	if (std::regex_search(f,sm,ex) && sm.size()>1){
	  //std::cout << sm.str(1) <<std::endl;
	  Index=stoul(sm.str(1).c_str());
	}
	else {
	  std::cout << "Couldn't extract index from filename " << f <<std::endl;
	  return 0;
	}
      }
      else {
	//just use order of entry
	Index++;
      }
      //now open unitcell file
      std::ifstream infile;
      infile.open(f.c_str(),ios::in | ios::binary);
      if (infile.is_open()){
	std::cout << "Loading " << f << std::endl;
	ajaj::UnitCell AA(ajaj::load_UnitCell_binary(infile,charge_rules,basis));
	std::complex<double> overlap(ajaj::Overlap(AA,AA));
	if (MElements.size() && basis.size()!=dim){
	  std::cout << "UnitCell Basis size doesn't match operator dimensions! " << basis.size() << " " << dim << std::endl;
	  return 0;
	}

	//measure.... //how many operators, separation etc...
	indexed_results.emplace_back(Index,ajaj::Data(abs(overlap)));

      }
      else {
	std::cout << "Couldn't open " << f << std::endl;
	return 0;
      }
    }

    //sort

    if (RuntimeArgs.use_filename_index()){
      //lambda expression
      std::sort(indexed_results.begin(),indexed_results.end(),[] (std::pair<size_t, ajaj::Data> a, std::pair<size_t, ajaj::Data> b) {
	  return b.first > a.first;
	});
    }

    ajaj::DataOutput results_file("UnitCellResults.dat","Index,abs(Overlap)");

    for (auto&& i : indexed_results){
      results_file.push(i.first,i.second);
    }


    return 1;
  }
  return 0;
}
