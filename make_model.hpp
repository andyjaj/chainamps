/** @file make_model.hpp
 * use args to direct correct model construction
 */
#ifndef MAKE_MODEL_H
#define MAKE_MODEL_H

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>

#include "command_line_input.hpp"
#include "common_defs.hpp"
#include "model.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "sparse_interface.hpp" //for entry() and SparseMatrix

//builtin models
#include "vertex_generators/continuum_Ising_generator.hpp"
//#include "vertex_generators/continuum_Ising_generator_no_sector.hpp"
#include "vertex_generators/continuum_ff_generator.hpp"
#include "vertex_generators/ll_generator.hpp"
#include "vertex_generators/spin_half_generator.hpp"
#include "vertex_generators/old_xxx_generator.hpp"

namespace ajaj{

  enum class BuiltinModels : std::size_t {
    cm_ising,
    cm_ff, 
    xxx,
    ll, 
    spin_half,
    user,
    last
  };

  struct NameGroup{
    BuiltinModels int_name;
    std::vector<std::string> name_variants;
    Vertex (*generator) (const VertexParameterArray&);
    MPO_matrix (*makeH) (const Vertex&, const VertexParameterArray&);

    NameGroup(BuiltinModels i, std::vector<std::string> n, Vertex (*g) (const VertexParameterArray&),MPO_matrix (*h) (const Vertex&, const VertexParameterArray&)) : int_name(i), name_variants(n),generator(g),makeH(h) {};
  };

  //recognised models, or user defined
  static const NameGroup m_names[7] = {
    {BuiltinModels::cm_ising, {"continuum_Ising", "continuum_Ising", "continuum ising", "Continuum_Ising", "Continuum Ising", "CONTINUUM_ISING", "CONTINUUM ISING"},&continuumIsing::VertexGenerator,&continuumIsing::MakeHamiltonian},
    {BuiltinModels::cm_ff, {"continuum_free_fermion", "continuum free fermion", "Continuum_Free_Fermion", "Continuum Free Fermion", "CONTINUUM_FREE_FERMION","CONTINUUM FREE FERMION"},&continuumff::VertexGenerator,&continuumff::MakeHamiltonian},
    {BuiltinModels::xxx, {"OLD_XXX"},&oldxxx::VertexGenerator,&oldxxx::MakeHamiltonian},
    {BuiltinModels::ll, {"ll","LL","Luttinger_liquid","luttinger_liquid","Luttinger liquid", "luttinger liquid", "LUTTINGER_LIQUID","LUTTINGER LIQUID"},&ll::VertexGenerator,&ll::MakeHamiltonian},
    {BuiltinModels::spin_half, {"spin half","SPIN_HALF","SPIN HALF"},&spin_half::VertexGenerator,&spin_half::MakeHamiltonian},
    {BuiltinModels::user, {"User_Defined","user_defined","User Defined","user defined","USER","user","USER_DEFINED", "USER DEFINED" },nullptr,nullptr},
    {BuiltinModels::last, {""},nullptr,nullptr}
  }; 

  NameGroup ModelNameChecker(const std::string& ns){
    for (std::size_t b=static_cast<std::size_t>(BuiltinModels::cm_ising);b<static_cast<std::size_t>(BuiltinModels::last);++b){
      for (auto&& m: m_names[b].name_variants){
	if (m==ns) return m_names[b];
      }
    }
    return m_names[static_cast<std::size_t>(BuiltinModels::last)];
  }

  struct OpBlockInfo{
    std::string row_block;
    std::string col_block;
    std::complex<double> value;
    OpBlockInfo(std::string r, std::string c, std::complex<double> v) : row_block(r), col_block(c), value(v) {}
  };

  Model MakeUserModel(const std::string& vertex_basis_filename,const std::string& vertex_hamiltonian_filename,const std::string& vertex_operators_filename,const CouplingArray& couplings);
  MPO_matrix MakeGeneralHMPO(const Vertex& v, const CouplingArray& couplings);

  Model MakeModelFromArgs(Base_Args& cmdln) {
    //check valid
    if (cmdln.is_valid()){
      std::ifstream infile;
      infile.open(cmdln.filename().c_str(),std::ios::in);
      if (infile.is_open()){
	std::cout << "Opened " << cmdln.filename() << std::endl;

	//read input file line by line into buffer
	std::string line;
	std::vector<std::string> stringbuffer;
	while (getline(infile,line)){
	  if (!line.empty())
	    stringbuffer.push_back(line);
	}
	infile.close();

	if (stringbuffer.size()!=3){
	  std::cout << "Input file error. File needs 3 non empty lines, but " << stringbuffer.size() << " present." << std::endl;
	}
	else { //parse strings
	  //first line is Model Name
	  const std::string& ModelName(stringbuffer.at(0));
	  NameGroup the_model(ModelNameChecker(ModelName));
	  VertexParameterArray vp;
	  //old style, needs killing off
	  VertexParameterArray cp;
	  CouplingArray cpa;

	  if (the_model.int_name==BuiltinModels::last){
	    std::cout << "No builtin model with name " << ModelName << std::endl;
	    //will default to empty model at end
	  }
	  else if(the_model.int_name==BuiltinModels::user) {
	    std::cout << "User Defined" <<std::endl;
	    //file format is model name
	    //spectrum and matrix element files (instead of vertex params)
	    std::istringstream iss1(stringbuffer.at(1));
	    std::string Vertex_Basis_Filename;
	    std::string Vertex_Hamiltonian_Filename;
	    std::string Vertex_Operators_Filename;
	    if (iss1 >> Vertex_Basis_Filename && iss1 >> Vertex_Hamiltonian_Filename && iss1 >> Vertex_Operators_Filename ) { //better be 3 filenames on this line
	      //also read in coupling params
	      std::istringstream iss2(stringbuffer.at(2));
	      Coupling c;
	      while (iss2 >> c){
		cpa.push_back(c);
	      }
	      //print error if no couplings, but proceed
	      if (cpa.size()==0){
		std::cout << "No inter vertex couplings defined!" << std::endl;
	      }
	      //what if couplings didn't read in correctly?
	      //won't be able to find them and will error below
	      return MakeUserModel(Vertex_Basis_Filename,Vertex_Hamiltonian_Filename,Vertex_Operators_Filename,cpa);
	    }
	    //didn't work
	    std::cout << "Not enough files specified." <<std::endl;
	    return Model();
	  }

	  else {
	    //builtin model
	    //now parse next line, containing vertex params
	    {
	      std::istringstream iss1(stringbuffer.at(1));
	      std::string word;
	      while (iss1 >> word){
		//use word
		std::string paramname(word);
		//advance
		if (!(iss1 >> word)) break;
		double paramvalue(stod(word));
		vp.push_back(VertexParameter(paramname,paramvalue));
		//vp.back().print();
	      }
	    }
	    //last line, containing coupling params	    
	    {
	      std::istringstream iss2(stringbuffer.at(2));
	      std::string word;
	      while (iss2 >> word){
		//use word
		std::string paramname(word);
		//advance
		if (!(iss2 >> word)) break;
		double paramvalue(stod(word));
		cp.push_back(VertexParameter(paramname,paramvalue));
	      }
	    }
	    std::cout << std::endl;
	    std::cout << "Builtin model, coupled array of: " << the_model.name_variants.back() << " vertices" << std::endl;
	    std::cout << "Vertex Parameters:" <<std::endl;
	    for (auto&& p : vp){
	      p.print();
	    }
	    std::cout << "Inter Vertex Coupling Parameters:" <<std::endl;
	    for (auto&& p : cp){
	      p.print();
	    }
	    if (cp.size()==0){
	      std::cout << "No inter vertex couplings defined! Are you sure you meant for this? Possible format error in input file.";
	    }
	    std::cout << std::endl;
	    

	    return Model(vp,the_model.generator,cp,the_model.makeH);
	  }
	}
      }
      else {
	std::cout << "Could not open " << cmdln.filename() << std::endl;
      }
    }
    std::cout << "Invalid command line arguments." << std::endl;
    std::cout << "Returning empty model..." << std::endl;
    return Model();
  }

  //user defined model, using file input
  //uses NRVO
  //really important not to mess this up, as it can invalidate the ref to ChargeRules used by all the states.
  Model MakeUserModel(const std::string& vertex_basis_filename,const std::string& vertex_hamiltonian_filename,const std::string& vertex_operators_filename,const CouplingArray& couplings){
    //construct model_vertex
    //read in basis file
    ajaj::Model UserModel;
    std::ifstream basisfile;
    basisfile.open(vertex_basis_filename.c_str(),std::ios::in);
    if (basisfile.is_open()){
      //process first row for charge defs
      //use getline, stringstream, stoi
      //ajaj::Vertex UserVertex;
      UserModel.vertex.Name=std::string("User Defined Vertex");
      std::string basisfileline;
      if (getline(basisfile,basisfileline)){
	std::istringstream basis_iss(basisfileline);
	std::string word;
	if (basis_iss >> word) {
	  if (word[0]=='D') { //starts with a D for definition line
	    int readinteger;
	    while (basis_iss >> readinteger) {
	      UserModel.vertex.ChargeRules.push_back(static_cast<QuantumNumberInt>(readinteger));
	    }
	  }
	  else if (word[0]=='#') {
	    //format is Z_n or Z
	    while (basis_iss >> word){
	      int tempint(word.find("_")!=std::string::npos ?  stoi(word.substr(word.rfind("_")+1,word.length())) : 0);
	      if (tempint< std::numeric_limits<short int>::max() && tempint > std::numeric_limits<short int>::min()){
		UserModel.vertex.ChargeRules.push_back(static_cast<QuantumNumberInt>(tempint));
	      }
	      else {
		std::cout << "Invalid charge definition!" <<std::endl;
		UserModel=Model(); 
		return UserModel;
	      }
	    }
	
	  }
	  else { //no definition
	    std::cout << "No definitions of n for Z_n charges!" <<std::endl;
	    UserModel=Model(); 
	    return UserModel;
	  }
	}
      }

      while (getline(basisfile,basisfileline)){
	std::istringstream basis_iss(basisfileline);
	int readinteger;
	if (basis_iss >> readinteger){ //dummy index
	  //read all charges
	  ajaj::QNVector charge_values;
	  while (basis_iss >> readinteger) {
	    charge_values.push_back(static_cast<QuantumNumberInt>(readinteger));
	  }
	  UserModel.vertex.Spectrum.push_back(ajaj::EigenState(UserModel.vertex.ChargeRules,charge_values,0.0));
	}
      }
      basisfile.close();

      std::cout << "MODEL'S LOCAL BASIS" <<std::endl;
      UserModel.basis().print();

      std::ifstream operatorsfile;
      operatorsfile.open(vertex_operators_filename.c_str(),std::ios::in);
      if (operatorsfile.is_open()){
	//read in operator names
	std::string operatorsfileline;
	if (getline(operatorsfile,operatorsfileline)){
	  std::istringstream operators_iss(operatorsfileline);
	  std::string word;
	  while (operators_iss >> word) {
	    UserModel.vertex.Operators.push_back(VertexOperator(word,UserModel.vertex.basis().size()));
	  }
	}
	//next M flag
	std::vector<bool>is_Measured;
	if (getline(operatorsfile,operatorsfileline)){
	  std::istringstream operators_iss(operatorsfileline);
	  std::string Mflag;
	  while (operators_iss >> Mflag) {
	    if (Mflag=="M"){
	      is_Measured.push_back(1);
	    }
	    else {
	      is_Measured.push_back(0);
	    }
	  }
	}
	if (is_Measured.size()!=UserModel.vertex.Operators.size()){std::cout << "Wrong number of Measured flags specified" << std::endl; UserModel=Model(); return UserModel;}

	//now entries
	while (getline(operatorsfile,operatorsfileline)){
	  std::istringstream operators_iss(operatorsfileline);
	  int row;
	  int col;
	  std::complex<double> value;
	  if (operators_iss >> row && operators_iss >> col){
	    size_t increment(0);
	    while (operators_iss >> value){
	      //std::cout << value << " " ;
	      if (value!=0.0){
		UserModel.vertex.Operators.at(increment).MatrixElements.entry(row,col,value);
	      }
	      ++increment;
	    }
	    if (increment!= UserModel.vertex.Operators.size()){
	      std::cout << "Incorrect number of operator values specified, " << increment << " : " << UserModel.vertex.Operators.size() << std::endl; UserModel=Model(); return UserModel;
	    }  
	  }
	}
	operatorsfile.close();
	//finalise operators
	for (auto&& O : UserModel.vertex.Operators){
	  O.MatrixElements.finalise();
	  //O.MatrixElements.print();
	}

	//open vertex_hamiltonian,
	std::ifstream hamiltonianfile;
	hamiltonianfile.open(vertex_hamiltonian_filename.c_str(),std::ios::in);
	if (hamiltonianfile.is_open()){
	  UserModel.vertex.Operators.push_back(VertexOperator("Vertex_Hamiltonian",UserModel.vertex.basis().size()));
	  is_Measured.push_back(0);
	  SparseMatrix& h_array=UserModel.vertex.Operators.back().MatrixElements;
	  //ajaj::SparseMatrix h_array(UserModel.vertex.Spectrum.size(),UserModel.vertex.Spectrum.size());
	  std::string hfileline;
	  while (getline(hamiltonianfile,hfileline)){
	    std::istringstream hss(hfileline);
	    MPXInt row;
	    MPXInt col;
	    std::complex<double> value;
	    if (hss >> row && hss >> col && hss >> value){
	      if (value!=0.0){
		h_array.entry(row,col,value);
	      }
	    }
	    else {
	      std::cout << "Malformed vertex hamiltonian" << std::endl;
	      UserModel=Model();
	      return UserModel;
	    }
	  }
	  h_array.finalise();
	  //build MPO
	  UserModel.H_MPO=MakeGeneralHMPO(UserModel.vertex,couplings);
	  std::cout << "MPO matrix info:" <<std::endl;
	  UserModel.H_MPO.print_indices();
	  return UserModel;

	}
	else {
	  std::cout << "Couldn't open vertex hamiltonian file" <<std::endl;
	  UserModel=Model(); 
	  return UserModel;
	}
      }
      else {
	std::cout << "Couldn't open operators file" <<std::endl;
	UserModel=Model();
	return UserModel;
      }
    }
    else {
      std::cout << "Couldn't open basis file" << std::endl;
    }
    return UserModel;
   }

  MPO_matrix MakeGeneralHMPO(const Vertex& v, const CouplingArray& couplings){
    //v has the basis info and the operators, including  the vertex hamiltonian
    //couplings tells us which to incude in MPO, and with what parameters

    //first we assemble row and col lists of used coupling operators
    std::vector<OpBlockInfo> CouplingBlocks;

    for (auto&& c : couplings){
      std::cout << "Coupling: " << c <<std::endl; 
      CouplingBlocks.emplace_back(c.Op1Name,c.Op2Name,c.Value);

    }
    QNCombinations differencecombinations(v.basis(),1); //1 means use difference

    //std::cout << "Diff pairs: " << differencecombinations.InvolutionPairs.size() << " " << differencecombinations.OrderedPairs.size() << std::endl;
    //differencecombinations.print();

    MPXInt lineardim((2+CouplingBlocks.size()*differencecombinations.size())*v.basis().size());

    SparseMatrix harray(lineardim,lineardim,(2+couplings.size())*v.basis().size());
    MPXInt offset_to_last_block((1+CouplingBlocks.size()*differencecombinations.size())*v.basis().size());
    //std::cout << "offset " << offset_to_last_block<<std::endl;
    //start with the really easy bits, the Identities
    for (size_t i=0;i<v.basis().size();++i){
      harray.entry(i,i,1.0);
      harray.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
    }

    //now vertex_hamiltonian
    const SparseMatrix& vh(v.get_operator_matrix("Vertex_Hamiltonian"));
    for (size_t col=0; col<vh.cols();++col){
      for (size_t p=vh.get_p(col); p<vh.get_p(col+1);++p){
	harray.entry(offset_to_last_block+vh.get_i(p),col,vh.get_x(p));
      }
    }

    //now the complicated bits....

    //loop through each operator pair

    for (auto b=0;b<CouplingBlocks.size();++b){
      const SparseMatrix& row_sp=v.get_operator_matrix(CouplingBlocks[b].row_block);
      const SparseMatrix& col_sp=v.get_operator_matrix(CouplingBlocks[b].col_block);
      const bool sameflag(CouplingBlocks[b].row_block==CouplingBlocks[b].col_block);
      const MPXInt row_block_row_offset=(1+b*differencecombinations.size())*v.basis().size();
      const MPXInt row_block_col_offset=0;
      const MPXInt col_block_row_offset=offset_to_last_block;
      const MPXInt col_block_col_offset=(1+b*differencecombinations.size())*v.basis().size();

      //do row array first
      for (size_t col=0; col<row_sp.cols();++col){
	for (size_t p=row_sp.get_p(col); p<row_sp.get_p(col+1);++p){
	  MPXInt i=row_sp.get_i(p);
	  State diffstate=v.basis()[i]-v.basis()[col];
	  MPXInt row_block_d_offset(0);
	  MPXInt col_block_d_offset(0);
	  if (diffstate==-diffstate){ //involution block
	    for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	      if (diffstate==differencecombinations.InvolutionPairs[l].PairState){
		row_block_d_offset=l*v.basis().size();
		if (sameflag)
		  col_block_d_offset=l*v.basis().size();
		break;
	      }
	    }
	  }
	  else {
	    for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	      if (diffstate==-differencecombinations.OrderedPairs[l].PairState){
		row_block_d_offset=(differencecombinations.InvolutionPairs.size()+l)*v.basis().size();
		if (sameflag)
		  col_block_d_offset=(differencecombinations.size()-l-1)*v.basis().size();		
		break;
	      }
	    }
	  }

	  harray.entry(i+row_block_row_offset+row_block_d_offset,col,row_sp.get_x(p));
	  if (sameflag)
	    harray.entry(i+col_block_row_offset,col+col_block_col_offset+col_block_d_offset,CouplingBlocks[b].value*row_sp.get_x(p));
	}
      }

      if (!sameflag) {//need to do col_block parts
	for (size_t col=0; col<col_sp.cols();++col){
	  for (size_t p=col_sp.get_p(col); p<col_sp.get_p(col+1);++p){
	    MPXInt i=col_sp.get_i(p);
	    State diffstate=v.basis()[i]-v.basis()[col];
	    MPXInt col_block_d_offset(0);
	    if (diffstate==-diffstate){ //involution block
	      for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
		if (diffstate==differencecombinations.InvolutionPairs[l].PairState){
		  col_block_d_offset=l*v.basis().size();
		  break;
		}
	      }
	    }
	    else {
	      for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
		if (diffstate==-differencecombinations.OrderedPairs[l].PairState){
		  col_block_d_offset=(differencecombinations.size()-l-1)*v.basis().size();		
		  break;
		}
	      }
	    }
	    harray.entry(i+col_block_row_offset,col+col_block_col_offset+col_block_d_offset,CouplingBlocks[b].value*col_sp.get_x(p));
	  }
	}
      }
    }

    //assemble charges
    StateArray s;
    s.push_back(State(v.basis().getChargeRules())); //identity is block diagonal in charges
    for (auto&& c : CouplingBlocks){
      for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	s.push_back(differencecombinations.InvolutionPairs[l].PairState);
      }
      for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	s.push_back(differencecombinations.OrderedPairs[l].PairState);
      }
    }
    s.push_back(State(v.basis().getChargeRules())); //vertex_hamiltonian must be block diagonal in charges
    
    std::vector<MPXIndex> indices;
    indices.push_back(MPXIndex(1,v.basis())); //sigma primed
    indices.push_back(MPXIndex(1,s)); //b_left
    indices.push_back(MPXIndex(0,v.basis())); //sigma
    indices.push_back(MPXIndex(0,s)); //b_right
    
    return MPO_matrix(v.basis(),indices,harray.finalise()); //finalise and MOVE array
  }
  

}

#endif
