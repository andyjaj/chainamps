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
#include "vertex_generators/continuum_ff_generator.hpp"
#include "vertex_generators/ll_generator.hpp"
#include "vertex_generators/llsc_generator.hpp"
#include "vertex_generators/old_xxx_generator.hpp"

namespace ajaj{

  enum class BuiltinModels : std::size_t {
    cm_ising,
      cm_ff, 
      xxx,
      ll,
      llsemi, 
      user,
      last
      };

  struct NameGroup{
    BuiltinModels int_name;
    std::vector<std::string> name_variants;
    Vertex (*generator) (const VertexParameterArray&);
    MPO_matrix (*makeH) (const Vertex&, const CouplingArray&);

    NameGroup(BuiltinModels i, std::vector<std::string> n, Vertex (*g) (const VertexParameterArray&),MPO_matrix (*h) (const Vertex&, const CouplingArray&)) : int_name(i), name_variants(n),generator(g),makeH(h) {};
  };

  //recognised models, or user defined
  static const NameGroup m_names[7] = {
    {BuiltinModels::cm_ising, {"continuum_Ising", "continuum_Ising", "continuum ising", "Continuum_Ising", "Continuum Ising", "CONTINUUM_ISING", "CONTINUUM ISING"},&continuumIsing::VertexGenerator,&continuumIsing::MakeHamiltonian},
    {BuiltinModels::cm_ff, {"continuum_free_fermion", "continuum free fermion", "Continuum_Free_Fermion", "Continuum Free Fermion", "CONTINUUM_FREE_FERMION","CONTINUUM FREE FERMION"},&continuumff::VertexGenerator,&continuumff::MakeHamiltonian},
    {BuiltinModels::xxx, {"OLD_XXX"},&oldxxx::VertexGenerator,&oldxxx::MakeHamiltonian},
    {BuiltinModels::ll, {"ll","LL","Luttinger_liquid","luttinger_liquid","Luttinger liquid", "luttinger liquid", "LUTTINGER_LIQUID","LUTTINGER LIQUID"},&ll::VertexGenerator,&ll::MakeHamiltonian},
    {BuiltinModels::llsemi, {"llsc","LLSC","Luttinger_liquid_semiclassical", "LUTTINGER_LIQUID_SEMICLASSICAL","LUTTINGER LIQUID SEMICLASSICAL"},&llsemi::VertexGenerator,&llsemi::MakeHamiltonian},
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

  Model MakeUserModel(const std::string& vertex_basis_filename,const std::string& vertex_hamiltonian_filename,const std::vector<std::string>& vertex_operator_filenames, const std::vector<CouplingArray>& couplings, const std::vector<double>& times);

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

	if (stringbuffer.size()<3){
	  std::cout << "Input file error. File needs 3 non empty lines, but " << stringbuffer.size() << " present." << std::endl;
	}
	else { //parse strings
	  if (stringbuffer.size()>3){
	    std::cout << "Input file has more than 3 non empty lines." << std::endl;
	    std::cout << "Only the first coupling definition line will be used by static routines." << std::endl;
	    std::cout << "Time routines will interpret suitably formatted extra lines as time dependent coupling information." << std::endl;
	  }
	  //first line is Model Name
	  const std::string& ModelName(stringbuffer.at(0));
	  NameGroup the_model(ModelNameChecker(ModelName));
	  VertexParameterArray vp;
	  
	  if (the_model.int_name==BuiltinModels::last){
	    std::cout << "No builtin model with name " << ModelName << std::endl;
	    //will default to empty model at end
	  }
	  else if(the_model.int_name==BuiltinModels::user) {
	    std::cout << "User Defined" <<std::endl;
	    //file format is model name
	    //spectrum and matrix element files (instead of vertex params)
	    std::istringstream iss(stringbuffer.at(1));
	    std::string Vertex_Basis_Filename;
	    std::string Vertex_Hamiltonian_Filename;
	    std::vector<std::string> Vertex_Operator_Filenames;
	    if (iss >> Vertex_Basis_Filename && iss >> Vertex_Hamiltonian_Filename){
	      std::string op_filename;
	      while (iss >> op_filename){
		Vertex_Operator_Filenames.push_back(op_filename);
	      }
	      //also read in coupling params
	      std::vector<CouplingArray> cpas;
	      std::vector<double> times;
	      //do we have time dep data, or not?
	      for (auto c_idx=2;c_idx<stringbuffer.size();++c_idx){
		std::istringstream css(stringbuffer.at(c_idx));
		css >> std::ws;
		double temp_time;
		bool is_time(0);
		if (std::isdigit(css.peek())){
		  is_time=1;
		  css >> temp_time;
		}
		else if (stringbuffer.size()>3){ //no time dep info and yet more than 1 coupling line
		  std::cout << "ERROR: no time data, but more than one coupling defs line!" << std::endl;
		  return Model();
		}
		cpas.push_back(CouplingArray());
		Coupling c;
		while (css >> c){
		  cpas.back().push_back(c);
		}
		if (cpas.back().size() && is_time){
		  times.push_back(temp_time);
		}
	      }
	      //print error if no couplings, but proceed
	      if (cpas.size()==0 || cpas[0].size()==0){
		std::cout << "NO INTER-VERTEX COUPLINGS DEFINED!" << std::endl;
		std::cout << "This is allowed, but is it what you intended?" <<std::endl;
	      }
	      //what if couplings didn't read in correctly?
	      //won't be able to find them and will error below
	      return MakeUserModel(Vertex_Basis_Filename,Vertex_Hamiltonian_Filename,Vertex_Operator_Filenames,cpas,times);
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
	    //last lines, containing coupling params

	    std::vector<CouplingArray> cpas;
	    std::vector<double> times;

	    for (auto c_idx=2;c_idx<stringbuffer.size();++c_idx){
	      std::istringstream css(stringbuffer.at(c_idx));
	      css >> std::ws;
	      double temp_time;
	      bool is_time(0);
	      if (std::isdigit(css.peek())){
		is_time=1;
		css >> temp_time;
	      }
	      else if (stringbuffer.size()>3){ //no time dep info and yet more than 1 coupling line
		std::cout << "ERROR: no time data, but more than one coupling defs line!" << std::endl;
		return Model();
	      }
	      cpas.push_back(CouplingArray());
	      //Coupling c;
	      std::string word;
	      std::complex<double> value;
	      while (css >> word >> value){
		cpas.back().push_back(Coupling(word,value));
	      }
	      if (cpas.back().size() && is_time){
		times.push_back(temp_time);
	      }
	    }
	    //print error if no couplings, but proceed
	    if (cpas.size()==0 || cpas[0].size()==0){
	      std::cout << "NO INTER-VERTEX COUPLINGS DEFINED!" << std::endl;
	      std::cout << "This is allowed, but is it what you intended?" <<std::endl;
	    }
	  
	    std::cout << std::endl;
	    std::cout << "Builtin model, coupled array of: " << the_model.name_variants.back() << " vertices" << std::endl;
	    std::cout << "Vertex Parameters:" <<std::endl;
	    for (auto&& p : vp){
	      p.print();
	    }
	    if (cpas.size()){
	      std::cout << "Inter Vertex Coupling Parameters:" <<std::endl;
	      for (auto&& c : cpas[0]){
		std::cout << c <<std::endl;
	      }
	    }
	    
	    else {
	      std::cout << "No inter vertex couplings defined! Are you sure you meant for this? Possible format error in input file.";
	    }
	    std::cout << std::endl;
	    

	    return Model(vp,the_model.generator,cpas,the_model.makeH,times);
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
  Model MakeUserModel(const std::string& vertex_basis_filename,const std::string& vertex_hamiltonian_filename,const std::vector<std::string>& vertex_operator_filenames,const std::vector<CouplingArray>& couplings,const std::vector<double>& times=std::vector<double>()){
    //construct model_vertex
    //read in basis file
    ajaj::Model UserModel(couplings,&MakeGeneralHMPO,times);
    UserModel.vertex.Name=std::string("User Defined Vertex");

    uMPXInt bsize=0;
    try { bsize=stoul(vertex_basis_filename); }
    catch (const std::invalid_argument& ia){
      //not a number
    }
    if (bsize>0){
      UserModel.vertex.ChargeRules.push_back(static_cast<QuantumNumberInt>(0));
      UserModel.vertex.Spectrum=ajaj::Basis(ajaj::StateArray(bsize,ajaj::State(UserModel.vertex.ChargeRules)),std::vector<double>(bsize,0.0));
    }
    else {
      std::ifstream basisfile;
      basisfile.open(vertex_basis_filename.c_str(),std::ios::in);
      if (basisfile.is_open()){
	//process first row for charge defs
	//use getline, stringstream, stoi
	//ajaj::Vertex UserVertex;
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
	      if (UserModel.vertex.ChargeRules.size()==0){
		std::cout <<"No charge rules, assuming no conserved Abelian charges..." <<std::endl;
		UserModel.vertex.ChargeRules.push_back(static_cast<QuantumNumberInt>(0));
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
	
	int dummyinteger(0);
	while (getline(basisfile,basisfileline)){
	  std::istringstream basis_iss(basisfileline);
	  if (basis_iss >> dummyinteger){ //first is a dummy index
	    //read all charges
	    ajaj::QNVector charge_values;
	    int readinteger;
	    while (basis_iss >> readinteger) {
	      charge_values.push_back(static_cast<QuantumNumberInt>(readinteger));
	    }
	    if (charge_values.size())
	      UserModel.vertex.Spectrum.push_back(ajaj::EigenState(UserModel.vertex.ChargeRules,charge_values,0.0));
	  }
	}
	basisfile.close();
      }
      else {
	std::cout << "Couldn't open basis file" << std::endl;
	UserModel=Model(); 
	return UserModel;
      }
      
      if (!UserModel.vertex.Spectrum.size()){ //still zero?
	std::cout << "No states defined? Returning empty model." <<std::endl;
	UserModel=Model(); 
	return UserModel;
      }
    }
  
    std::cout << "MODEL'S LOCAL BASIS" <<std::endl;
    UserModel.basis().print();
    
    for (auto&& f : vertex_operator_filenames){
      
      std::cout <<"Opening operator file " << f <<std::endl;
      
      std::ifstream operatorfile;
      operatorfile.open(f.c_str(),std::ios::in);
      if (operatorfile.is_open()){
	size_t previous_num_ops(UserModel.vertex.Operators.size());
	//read in operator names
	std::string operatorfileline;
	if (getline(operatorfile,operatorfileline)){
	  std::istringstream operators_iss(operatorfileline);
	  std::string word;
	  while (operators_iss >> word) {
	    UserModel.vertex.Operators.push_back(VertexOperator(word,UserModel.vertex.basis().size()));
	  }
	}
	size_t ops_in_file(UserModel.vertex.Operators.size()-previous_num_ops);
	//next, possible dummy line
	if (!std::isdigit(operatorfile.peek())){
	  //std::vector<std::string>dummy_names;;
	  if (getline(operatorfile,operatorfileline)){
	    /*std::istringstream operators_iss(operatorfileline);
	      std::string name;
	      while (operators_iss >> name) {
	      dummy_names.push_back(name);
	      std::cout << name << std::endl;
	      }*/
	  }
	}

	//now entries
	while (getline(operatorfile,operatorfileline)){
	  std::istringstream operators_iss(operatorfileline);
	  int row;
	  int col;
	  std::complex<double> value;
	  if (operators_iss >> row && operators_iss >> col){
	    size_t increment(0);
	    while (operators_iss >> value){
	      //std::cout << value << " " ;
	      if (value!=0.0){
		UserModel.vertex.Operators.at(increment+previous_num_ops).MatrixElements.entry(row,col,value);
	      }
	      ++increment;
	    }
	    if (increment!= ops_in_file){
	      std::cout << "Incorrect number of operator values specified, " << increment << " : " << UserModel.vertex.Operators.size() << std::endl; UserModel=Model(); return UserModel;
	    }
	  }
	}
	operatorfile.close();
      }
      else {
	std::cout << "Couldn't open operators file" <<std::endl;
	UserModel=Model();
	return UserModel;
      }
    }
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
      SparseMatrix& h_array=UserModel.vertex.Operators.back().MatrixElements;
      //ajaj::SparseMatrix h_array(UserModel.vertex.Spectrum.size(),UserModel.vertex.Spectrum.size());
      std::string hfileline;
      while (getline(hamiltonianfile,hfileline)){
	std::istringstream hss(hfileline);
	MPXInt row;
	MPXInt col;
	std::complex<double> value;
	hss >> std::ws;
	if (hss.peek()==EOF){
	  continue;
	}
	else if (hss >> row && hss >> col && hss >> value){
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
      UserModel.H_MPO=MakeGeneralHMPO(UserModel.vertex,couplings.at(0));
      std::cout << "MPO matrix info:" <<std::endl;
      UserModel.H_MPO.print_indices();
      UserModel.H_MPO.print_sparse_info();
      return UserModel;

    }
    else {
      std::cout << "Couldn't open vertex hamiltonian file" <<std::endl;
      UserModel=Model(); 
      return UserModel;
    }
    return UserModel;
  }

  MPO_matrix MakeGeneralHMPO(const Vertex& v, const CouplingArray& couplings){
    //v has the basis info and the operators, including  the vertex hamiltonian
    //couplings tells us which to incude in MPO, and with what parameters

    CouplingArray InterVertexCouplings;
    CouplingArray LocalTerms;

    for (auto&& c : couplings){
      if (c.two_vertex()){
	std::cout << "Inter vertex coupling: " << c <<std::endl; 
	InterVertexCouplings.push_back(c);
      }
      else {
	std::cout << "Single vertex term: " << c <<std::endl; 
	LocalTerms.push_back(c);
      }
    }

    QNCombinations differencecombinations(v.basis(),1); //1 means use difference

    MPXInt lineardim=(2+InterVertexCouplings.size()*differencecombinations.size())*v.basis().size();
    SparseMatrix harray(lineardim,lineardim,(2+InterVertexCouplings.size())*v.basis().size());
    MPXInt offset_to_last_block((1+InterVertexCouplings.size()*differencecombinations.size())*v.basis().size());
    //std::cout << "offset " << offset_to_last_block<<std::endl;
    //start with the really easy bits, the Identities
    for (size_t i=0;i<v.basis().size();++i){
      harray.entry(i,i,1.0);
      harray.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
    }

    //now vertex_hamiltonian
    //plus any local terms...
    {
      const SparseMatrix& vh(v.get_operator_matrix("Vertex_Hamiltonian"));
      for (size_t col=0; col<vh.cols();++col){
	for (size_t p=vh.get_p(col); p<vh.get_p(col+1);++p){
	  harray.entry(offset_to_last_block+vh.get_i(p),col,vh.get_x(p));
	}
      }
    }

    for (auto&& l : LocalTerms){
      const SparseMatrix& op(v.get_operator_matrix(l.Op1Name));
      for (size_t col=0; col<op.cols();++col){
	for (size_t p=op.get_p(col); p<op.get_p(col+1);++p){
	  harray.entry(offset_to_last_block+op.get_i(p),col,l.Value*op.get_x(p));
	}
      }
    }

    //now the complicated bits....

    //loop through each operator pair

    for (auto b=0;b<InterVertexCouplings.size();++b){
      //const SparseMatrix& row_sp=v.get_operator_matrix(CouplingBlocks[b].row_block);
      //const SparseMatrix& col_sp=v.get_operator_matrix(CouplingBlocks[b].col_block);
      //const bool sameflag(CouplingBlocks[b].row_block==CouplingBlocks[b].col_block);

      const SparseMatrix& row_sp=v.get_operator_matrix(InterVertexCouplings[b].Op1Name);
      const SparseMatrix& col_sp=v.get_operator_matrix(InterVertexCouplings[b].Op2Name);
      const bool sameflag(InterVertexCouplings[b].Op1Name==InterVertexCouplings[b].Op2Name);

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
	    harray.entry(i+col_block_row_offset,col+col_block_col_offset+col_block_d_offset,InterVertexCouplings[b].Value*row_sp.get_x(p));
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
	    harray.entry(i+col_block_row_offset,col+col_block_col_offset+col_block_d_offset,InterVertexCouplings[b].Value*col_sp.get_x(p));
	  }
	}
      }
    }

    //assemble charges
    StateArray s;
    s.push_back(State(v.basis().getChargeRules())); //identity is block diagonal in charges
    for (auto&& c : InterVertexCouplings){
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
