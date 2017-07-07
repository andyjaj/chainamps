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

#include "common_defs.hpp"
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

    //if we want to add time info, this is where we should do it
    std::vector<double> times;
    if (!RuntimeArgs.time_filename().empty()){
      std::ifstream time_file;
      time_file.open(RuntimeArgs.time_filename().c_str(),ios::in);
      
      if (time_file.is_open()){
	size_t idx;
	double time;
	std::string sdump;
	time_file >> std::ws;
	if (time_file.peek()=='#'){getline(time_file,sdump);}
	while (time_file >> idx >> time){
	  getline(time_file,sdump);
	  times.push_back(time);
	}
	
      }
      else {
	std::cout <<"Could not open specified time data file '" << RuntimeArgs.time_filename() << "'." <<std::endl;
	return 1;
      }
    }

    ajaj::Vertex iMEAS_vertex; // a dummy vertex
    std::vector<ajaj::ShiftedOperatorInfo> opinfo;

    //read in (possible) operators    

    if (RuntimeArgs.operator_filenames().size()) {
      for (auto&& opfn : RuntimeArgs.operator_filenames()){
	
	ajaj::ShiftedOperatorInfo this_opinfo(opfn);
	std::ifstream is;
	is.open(this_opinfo.Name.c_str(),std::ios::in);
	if (is.is_open()){
	  std::istringstream opss(this_opinfo.Name);
	  std::string opname;
	  getline(opss,opname,'.');
	  if (!iMEAS_vertex.operator_exists(opname)){
	    iMEAS_vertex.Operators.emplace_back(opname);
	  
	    iMEAS_vertex.Operators.back().MatrixElements=ajaj::load_SparseMatrix(is);
	    if (iMEAS_vertex.Operators.back().MatrixElements.rows()!=iMEAS_vertex.Operators.back().MatrixElements.cols()){
	      std::cout << "Malformed operator in " << this_opinfo.Name << std::endl;
	      return 0;
	    }
	  }
	  this_opinfo.Name=opname;
	  opinfo.push_back(this_opinfo);
	}
	else {
	  std::cout << "Couldn't open " << this_opinfo.Name << std::endl;
	  return 0;
	}
      }
    }

    if (RuntimeArgs.nev()) std::cout << "Number of transfer matrix eigenvalues requested: " << RuntimeArgs.nev() << std::endl;

    if (opinfo.size()){
      std::cout << "Measuring " << (RuntimeArgs.two_point() ? 2 : 1) << "-point function" << std::endl;
      std::cout << opinfo.front();
      if (RuntimeArgs.two_point()){
	std::cout << "(i) " << opinfo.back() << "(i+" << RuntimeArgs.separation() <<")";
      }
      std::cout <<std::endl;
    }
    //check dims
    ajaj::MPXInt dim(iMEAS_vertex.Operators.size() ? iMEAS_vertex.Operators[0].MatrixElements.rows() : 0);

    if (iMEAS_vertex.Operators.size()>1){
      for (auto&& m : iMEAS_vertex.Operators){
	if (dim!= m.MatrixElements.rows()){
	  std::cout << "Operator dimensions don't match! " << dim << " " << m.MatrixElements.rows() <<std::endl;
	  return 0;
	}
      }
    }

    std::vector<std::pair<size_t,ajaj::Data> > indexed_results;
    size_t Index(0);
    std::regex ex("_([0-9]+)\\.UNITCELL$"); //regex to match number before filename only.

    std::vector<ajaj::NamedMPO_matrix> OperatorMPOs; //if defined
    
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

      //are we using a time file, and does the unitcell index exceed its size?
      if (times.size() && Index<=times.size()){
      
      //now open unitcell file
      std::ifstream infile;
      infile.open(f.c_str(),ios::in | ios::binary);
      if (infile.is_open()){
	std::cout << "Loading " << f << std::endl;
	ajaj::UnitCell AA(ajaj::load_UnitCell_binary(infile,iMEAS_vertex.ChargeRules,iMEAS_vertex.Spectrum));//populates basis
	ajaj::State TargetState(iMEAS_vertex.ChargeRules,RuntimeArgs.target());

	for (auto&& m : AA.Matrices){
	  m.print_sparse_info();
	}
	indexed_results.emplace_back(Index,ajaj::Data());
	if (times.size()){
	  indexed_results.back().second.Real_measurements.emplace_back(times.at(Index-1));
	}
	if (RuntimeArgs.calc_entanglement()){
	  indexed_results.back().second.Real_measurements.emplace_back(AA.Entropy());
	}
	if (RuntimeArgs.vert_entanglement()){
	  indexed_results.back().second.Real_measurements.emplace_back(MultiVertexEntropy(AA,RuntimeArgs.vert_entanglement()));
	}
	if (RuntimeArgs.nev()){
	  std::vector<std::complex<double> > Transfer_eigs(ajaj::TransferMatrixEigs(AA,RuntimeArgs.nev(),TargetState));
	  for (auto&& eig : Transfer_eigs) {
	    indexed_results.back().second.Real_measurements.emplace_back(abs(eig));
	  }
	}
	if (iMEAS_vertex.Operators.size()){
	  if (iMEAS_vertex.Spectrum.size()!=dim){
	    std::cout << "UnitCell Basis size doesn't match operator dimensions! " << iMEAS_vertex.Spectrum.size() << " " << dim << std::endl;
	    return 0;
	  }
	  //measure.... //how many operators, separation etc...
	  if (!OperatorMPOs.size()){ //if haven't populated Operators yet, do it now
	    for (auto&& shiftedop : opinfo){
	      //make name
	      std::ostringstream mponame;
	      mponame << shiftedop.Name;
	      if (shiftedop.Factor!=0.0){
		mponame << "@" << shiftedop.WhichCharge << "@" << shiftedop.Factor;
 	      }
	      OperatorMPOs.emplace_back(mponame.str(),iMEAS_vertex.make_one_site_operator(shiftedop));
	    }
	  }

	  if (OperatorMPOs.size()==1){
	    if (RuntimeArgs.two_point()==0) //single site
	      indexed_results.back().second.Complex_measurements.emplace_back(OneVertexMeasurement(OperatorMPOs[0].Matrix,AA));
	    else //same operator, with separation
	      indexed_results.back().second.Complex_measurements.emplace_back(TwoVertexMeasurement(OperatorMPOs[0].Matrix,OperatorMPOs[0].Matrix,AA,RuntimeArgs.separation()));
	  }
	  else if (OperatorMPOs.size()==2){ //two (possibly different) operators
	    indexed_results.back().second.Complex_measurements.emplace_back(TwoVertexMeasurement(OperatorMPOs[0].Matrix,OperatorMPOs[1].Matrix,AA,RuntimeArgs.separation()));
	  }
	}
      }
      else {
	std::cout << "Couldn't open " << f << std::endl;
	return 0;
      }
      }
    }
    //sort

    if (RuntimeArgs.use_filename_index()){
      //lambda expression
      std::sort(indexed_results.begin(),indexed_results.end(),[] (std::pair<size_t, ajaj::Data> a, std::pair<size_t, ajaj::Data> b) {
	  return b.first > a.first;
	});
    }

    std::ostringstream mss;
    if (opinfo.size()){
      mss << opinfo.front();
      if (RuntimeArgs.two_point())
	mss << "_" << opinfo.back() << "_" << RuntimeArgs.separation();
    }

    std::ostringstream outfilename;
    outfilename << "UNITCELL_Results";
    if (!mss.str().empty())
      outfilename << "_" << mss.str();
    if (RuntimeArgs.nev()){
      outfilename << "_TransferMatrix_EVals_TargetState";
      for (auto t : RuntimeArgs.target()){
	outfilename << "_" << t;
      }
    }
    outfilename << ".dat";

    std::ostringstream commentstream;

    if (RuntimeArgs.nev()){
      commentstream << "Transfer Matrix EigenValues for Target State: ";
      for (auto t : RuntimeArgs.target()){
	commentstream << t <<" ";
      }
      commentstream << "\n";
    }

    commentstream << "Index";
    if (!RuntimeArgs.time_filename().empty()){
      commentstream <<",time";
    }
    if (RuntimeArgs.calc_entanglement()){
      commentstream << ",S_E";
    }
    if (RuntimeArgs.vert_entanglement()){
      commentstream << "," << RuntimeArgs.vert_entanglement() <<" vertex S_E";
    }

    for (size_t l=0; l<RuntimeArgs.nev();++l){
      commentstream << ",abs(Lambda_" << l+1 << ")";
    }
    if (opinfo.size()) {
      std::ostringstream opss;
      opss << "," << opinfo.front();
      if (RuntimeArgs.two_point()){
	opss << "(i)," << opinfo.back() << "(i+" << RuntimeArgs.separation() << ")";
      }
      commentstream << ",Re(" << opss.str() <<"),Im(" << opss.str() << ")";
    }

    ajaj::DataOutput results_file(outfilename.str(),commentstream.str());

    for (auto&& i : indexed_results){
      results_file.push(i.first,i.second);
    }

    std::cout << std::endl << "Done" << std::endl << "Results written to " << outfilename.str() <<std::endl;

    return 1;
  }
  return 0;
}
