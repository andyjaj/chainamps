/** @file command_line_input.hpp
 * Initial input from command line can be given in the form of a spectrum and matrix elements
 * A model still needs to be defined though for the Hamiltonian to be formed correctly.
 */
#ifndef COMMAND_LINE_INPUT_H
#define COMMAND_LINE_INPUT_H
#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include <vector>
#include <complex>
#include "optionparser/optionparser.h" //The Lean Mean C++ Option Parser. See .h file for license info.

namespace ajaj {

  std::istream& operator>>(std::istream& is, std::vector<short>& v)
  {
    //drop initial whitespace
    is >> std::ws;
    std::vector<short> temp;
    bool failure(0);
    std::string s;
    int n;

    while (getline(is,s,',')){ //reads to comma, or until eof

      n=stoi(s);

      if (n< std::numeric_limits<short int>::max() && n > std::numeric_limits<short int>::min()){
	temp.push_back(n);
      }
      else {
	failure=1;
	break;
      }
      if (is.eof()) break;
    }

    //getline fail?

    //std::cout << is.bad() << " " << is.fail() << " " << is.eof() << " " <<is.good() <<std::endl;

    if( failure )
      is.setstate(std::ios::failbit);
    else
      v=temp;

    return is;
  }

  typedef std::vector<std::pair<std::string,unsigned long> > StringIndexPairs;

  std::istream& operator>>(std::istream& is, StringIndexPairs& l)
  {
    //drop initial whitespace
    is >> std::ws;
    StringIndexPairs temp;
    bool failure(0);
    std::string s;

    while (getline(is,s,',')){ //reads to comma, or until eof
      //and get the next one too
      unsigned long n(0);
      std::string num;
      if (getline(is,num,',')){
	size_t i(0);
	n=stoul(num);
	temp.emplace_back(s,n);
      }
      else {
	failure=1;
	break;
      }
      if (is.eof()) break;
    }
    if( failure )
      is.setstate(std::ios::failbit);
    else
      l=temp;

    return is;
  }

  struct Arg: public option::Arg
  {
    static option::ArgStatus Unknown(const option::Option& option, bool msg)
    {
      if (msg) std::cout << "Unknown option '" << std::string(option.name,option.namelen) << "'" <<std::endl;
      return option::ARG_ILLEGAL;
    }

    static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
    {
      if (option.arg != 0 && option.arg[0] != 0 && option.arg[0] != '-')
	return option::ARG_OK;

      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a non-empty argument" <<std::endl;
      return option::ARG_ILLEGAL;
    }

    static option::ArgStatus PositiveNumeric(const option::Option& option, bool msg)
    {
      char* endptr = nullptr;
      if (option.arg != 0 && strtol(option.arg, &endptr, 10)>=0){};
      if (endptr != option.arg && *endptr == 0)
	return option::ARG_OK;

      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a positive numeric argument" <<std::endl;
      return option::ARG_ILLEGAL;
    }

    static option::ArgStatus PositiveDefiniteNumeric(const option::Option& option, bool msg)
    {
      char* endptr = nullptr;
      if (option.arg != 0 && strtol(option.arg, &endptr, 10)>0){};
      if (endptr != option.arg && *endptr == 0)
	return option::ARG_OK;

      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a positive definite numeric argument" <<std::endl;
      return option::ARG_ILLEGAL;
    }

    static option::ArgStatus ShortNumeric(const option::Option& option, bool msg)
    {
      char* endptr = nullptr;
      long n;
      if (option.arg != 0) n=strtol(option.arg, &endptr, 10); //if not null, then convert (base 10)
      if (endptr != option.arg && *endptr == 0 && n< std::numeric_limits<short int>::max() && n > std::numeric_limits<short int>::min())
	return option::ARG_OK;

      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a short integer as an argument" <<std::endl;
      return option::ARG_ILLEGAL;
    }
    static option::ArgStatus CommaSepShorts(const option::Option& option, bool msg)
    {
      if (option.arg != 0 && option.arg[0]){
	std::vector<short> testvec;
	std::istringstream ss(option.arg);
	if (ss >> testvec)
	  return option::ARG_OK;
      }  

      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a comma separated list of short integers as an argument" <<std::endl;
      return option::ARG_ILLEGAL;
    }
    static option::ArgStatus PositiveDouble(const option::Option& option, bool msg)
    {
      size_t pos(0);
      if (option.arg != nullptr && stod(option.arg,&pos)>0.0){
	  return option::ARG_OK;
      }
      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a positive numeric argument" <<std::endl;
      return option::ARG_ILLEGAL;
    }
    static option::ArgStatus FiniteMeasurementInfo(const option::Option& option, bool msg)
    {
      if (option.arg != 0 && option.arg[0] && option.arg[0] != '-'){
	StringIndexPairs check;
	std::istringstream ss(option.arg);
	ss >>check;
	if (check.size()) return option::ARG_OK;
      }

      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a comma separated list of measurements (operator,position,operator,position,...)" <<std::endl;
      return option::ARG_ILLEGAL;
    }

  };

  enum optionIndex {UNKNOWN,CHI,TRUNC,NUMBER_OF_STEPS,MINS,NUMBER_OF_EXCITED,NUMBER_OF_SWEEPS,WEIGHT_FACTOR,TROTTER_ORDER,TIME_STEPS,STEP_SIZE,MEASUREMENT_INTERVAL,INITIAL_STATE_NAME,SEPARATION,NOINDEX,OPERATORFILE,TARGET,FINITE_MEASUREMENT,NEV,ENTANGLEMENT,VERTEX_ENTANGLEMENT,C_SPECIFIER,TIMEFILE,FDMRG_MODE,IENERGY,SAVE_ALL};

  const option::Descriptor store_usage[2] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: STORE_OPERATORS.bin <model_filename>"
       "\tStores all the defined vertex operators for a model as csc SPARSEMATRIX files."},
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor iDMRG_usage[6] =
    {
      {UNKNOWN, 0,"", "",Arg::Unknown, "USAGE: iDMRG_DRV.bin [OPTIONS] <model_filename>"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, --bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
      {TRUNC,0,"e","truncation-error",Arg::PositiveDouble,"  -e <number>, --truncation-error=<number>"
       "  \tThe allowed truncation error, >= 0." },
      {NUMBER_OF_STEPS,0,"N","number-of-steps",Arg::PositiveNumeric,"  -N <number>, --number-of-steps=<number>"
       "  \tThe number of infinite volume steps, >= 0" },
      {TARGET,0,"T","target-charges",Arg::CommaSepShorts,"  -T <number>,<number>,<number>, --target-charges=<n>,<n>,<n>"
       "  \t Charges to target for the unit cell."},
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor fDMRG_usage[8] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: fDMRG_DRV.bin [OPTIONS] <model_filename> <number of vertices/chains> \n  <number of vertices/chains> must be EVEN.\n"
	"Options:"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, --bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
      {TRUNC,0,"e","truncation-error",Arg::PositiveDouble,"  -e <number>, --truncation-error=<number>"
       "  \tThe allowed truncation error >= 0." },
      {NUMBER_OF_EXCITED,0,"X","excited-states",Arg::PositiveNumeric,"  -X <number>, --excited-states=<number>"
       "  \tThe number of excited states to calculate, >= 0. If this is specified but a projective weight factor is not, then a default value of 100.0 is used." },
      {NUMBER_OF_SWEEPS,0,"F","finite-size-sweeps",Arg::PositiveNumeric,"  -F <number>, --finite-size-sweeps=<number>"
       "  \tThe number of finite size sweeps to perform, >= 0" },
      {WEIGHT_FACTOR,0,"W","weight-factor",Arg::PositiveDouble,"  -W <number>, --weight-factor=<number>"
       "  \tThe weight factor if calculating excited states. Must be > 0. Specifying this indicates that the number of requested excited states is at least 1." },
      {TARGET,0,"T","target-charges",Arg::CommaSepShorts,"  -T <number>,<number>,<number>, --target-charges=<n>,<n>,<n>"
       "  \tCharges for target state."},
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor iTEBD_usage[10] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: iTEBD_DRV.bin [-B <number> -n <number> -s <number> -O <number>] <model_filename>"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, \t--bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
      {TRUNC,0,"e","truncation-error",Arg::PositiveDouble,"  -e <number>, \t--truncation-error=<number>"
       "  \tThe allowed truncation error, >= 0." },
      {NUMBER_OF_STEPS,0,"n","time-steps",Arg::PositiveNumeric,"  -n <number>, \t--time-steps=<number>"
       "  \tThe number of time steps. Default is 1." },
      {STEP_SIZE,0,"s","step-size",Arg::PositiveDouble,"  -s <number>, \t--step-size=<number>"
       "  \tThe step size. Default is 0.1" },
      {TROTTER_ORDER,0,"O","trotter-order",Arg::PositiveNumeric,"  -O <number>, \t--trotter-order=<number>"
       "  \tThe Trotter order (currently 1 or 2). Second order (2) is default." },
      {MEASUREMENT_INTERVAL,0,"m","measurement-interval",Arg::PositiveNumeric,"  -m <number>, \t--measurement-interval=<number>"
       "  \tMeasurement at every <number> steps. Default is 1 (measurement at every step)."},
      {INITIAL_STATE_NAME,0,"i","initial-unit-cell",Arg::NonEmpty,"  -i <initial_unit_cell>, \t--initial-unit-cell=<initial_unit_cell>"
       "  \tSpecify an initial unit cell." },
      {C_SPECIFIER,0,"c","c-number-file",Arg::NonEmpty," -c <c-specifier-file>,"
       "  \tFile with c-numbers for initial state."},
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor TEBD_usage[12] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: TEBD_DRV.bin [OPTIONS] <model_filename> <number of vertices(chains)> \n  <number of vertices/chains> must be EVEN.\n"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, \t--bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
      {TRUNC,0,"e","truncation-error",Arg::PositiveDouble,"  -e <number>, \t--truncation-error=<number>"
       "  \tThe allowed truncation error, >= 0." },
      {NUMBER_OF_STEPS,0,"n","time-steps",Arg::PositiveNumeric,"  -n <number>, \t--time-steps=<number>"
       "  \tThe number of time steps. Default is 1." },
      {STEP_SIZE,0,"s","step-size",Arg::PositiveDouble,"  -s <number>, \t--step-size=<number>"
       "  \tThe step size. Default is 0.1" },
      {TROTTER_ORDER,0,"O","trotter-order",Arg::PositiveNumeric,"  -O <number>, \t--trotter-order=<number>"
       "  \tThe Trotter order (currently 1 or 2). Second order (2) is default." },
      {MEASUREMENT_INTERVAL,0,"m","measurement-interval",Arg::PositiveNumeric,"  -m <number>, \t--measurement-interval=<number>"
       "  \tMeasurement at every <number> steps. Default is 1 (measurement at every step)."},
      {INITIAL_STATE_NAME,0,"i","initial-state-name",Arg::NonEmpty,"  -i <initial_state_name>, \t--initial-state-name=<initial_state_name>"
       "  \tSpecify an initial state." },
      {FINITE_MEASUREMENT,0,"M","finite-measurement",Arg::FiniteMeasurementInfo,"  -M <opfile1>,<vertex1>[,<opfile2>,<vertex2>], \t--finite-measurement=<opfile1>,<vertex1>[,<opfile2>,<vertex2>]"
       "  \tSpecify a one or two point measurement."},
      {C_SPECIFIER,0,"c","c-number-file",Arg::NonEmpty,"  -c <c-specifier-file>,  \t--c-number-file=<c-specifier-file>"
       "  \tFile with c-numbers for initial state."},
      {SAVE_ALL,0,"A","save-all",Arg::None,"  -A, \t--save-all"
       "  \tSave files for all times. Also sets measurement-interval=1."},
      { 0, 0, 0, 0, 0, 0 }
    };  

  const option::Descriptor iMEAS_usage[11] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: UNITCELL_MEASURE.bin [OPTIONS] <unitcell_filename1> ... \n"},
      {OPERATORFILE,0,"O","operator filename",Arg::NonEmpty,"  -O <filename>, \t--operator-file=<filename>"
       "  \tFile containing sparse matrix definition of operator." },
      {SEPARATION,0,"S","separation",Arg::PositiveDefiniteNumeric,"  -S <number>, \t--separation=<number>"
       "\tDistance (in integer units) between two vertex operator measurements: Op1(0),Op2(<number>). Specifying any separation, including 0, is interpreted as a two point measurement (possibly squaring an operator)." },
      {NOINDEX,0,"X","No index",Arg::None,"  -X, \t--no-index"
       "  \tDon't extract index from filenames." },
      {NEV,0,"N","Number of eigenvalues",Arg::PositiveDefiniteNumeric,"  -N, \t--number-of-eigenvalues=<number>"
       "  \tNumber of transfer matrix eigenvalues to calculate. Requires a target to be set, using -T" },
      {TARGET,0,"T","target-charges",Arg::CommaSepShorts,"  -T <number>,<number>,<number>, \t--target-charges=<n>,<n>,<n>"
       "  \tCharges to target eigenvalues."},
      {ENTANGLEMENT,0,"E","Entanglement Entropy",Arg::None,"  -E, \t--entanglement-entropy"
       "  \tCalculate the bipartite entanglement for the infinite system."},
      {VERTEX_ENTANGLEMENT,0,"V","Multi Vertex Entanglement",Arg::PositiveDefiniteNumeric,"  -V <number>, \t--vertex-entanglement=<number>"
       "  \tCalculate the entanglement entropy for <number> consecutive vertices in the infinite system."},
      {TIMEFILE,0,"t","time-filename",Arg::NonEmpty," -t <filename>, \t--time-filename=<filename>"
       "  \tUse time data from <filename> to include times in output file."},
      {IENERGY,0,"H","Energy per Vertex",Arg::NonEmpty," -H <model_filename>, \t--hamiltonian-modelfile=<model_filename>"
       "  \tCalculate the energy per vertex for the unitcell using Hamiltonian defined in the model file."},
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor fMEAS_usage[4] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: FINITE_MEASURE.bin [OPTIONS] <model_filename> <number of vertices(chains)> <state_name1> ... \n  <number of vertices/chains> must be EVEN.\n"},
      {FINITE_MEASUREMENT,0,"M","finite-measurement",Arg::FiniteMeasurementInfo,"  -M <opfile1>,<vertex1>[,<opfile2>,<vertex2>], \t--finite-measurement=<opfile1>,<vertex1>[,<opfile2>,<vertex2>]"
       "  \tSpecify a one or two point measurement."},
      {FDMRG_MODE,0,"D","fDMRG-mode",Arg::None,"  -D, \t--fDMRG-mode"
       "  \tSpecial mode for fDMRG output files, needs no input filenames."},
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor TwoVE_usage[7] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: 2VE_DRV.bin [OPTIONS] <model_filename> \n Exact diagonalisation used to time evolve two vertices."},
      {STEP_SIZE,0,"s","step-size",Arg::PositiveDouble,"  -s <number>, \t--step-size=<number>"
       "  \tMeasurements are taken at time intervals separated by the step size. Default is 0.1" },
      {NUMBER_OF_STEPS,0,"n","time-steps",Arg::PositiveNumeric,"  -n <number>, \t--time-steps=<number>"
       "  \tThe number of time steps. Default is 1." },
      {INITIAL_STATE_NAME,0,"i","initial-state-name",Arg::NonEmpty,"  -i <initial_state_name>, \t--initial-state-name=<initial_state_name>"
       "  \tSpecify an initial state." },
      {FINITE_MEASUREMENT,0,"M","finite-measurement",Arg::FiniteMeasurementInfo,"  -M <opfile1>,<vertex1>[,<opfile2>,<vertex2>], \t--finite-measurement=<opfile1>,<vertex1>[,<opfile2>,<vertex2>]"
       "  \tSpecify a one or two point measurement."},
      {C_SPECIFIER,0,"c","c-number-file",Arg::NonEmpty," -c <c-specifier-file>,"
       "  \tFile with c-numbers for initial state."},
      { 0, 0, 0, 0, 0, 0 }
    };
  
  class Base_Args{

  protected:
    int argc_;
    char** argv_;
    const option::Descriptor* usage_;
    option::Stats stats;
    std::vector<option::Option> options;
    std::vector<option::Option> buffer;
    option::Parser parse;
    bool valid_;

  public:
    Base_Args(int argc, char* argv[], const option::Descriptor* usage) : argc_(argc-(argc>0)),argv_(argv+(argc>0)),usage_(usage),stats(usage_, argc_, argv_),options(stats.options_max),buffer(stats.buffer_max),parse(usage_, argc_, argv_, &options[0], &buffer[0]),valid_(0){
      std::cout<< "ChainAMPS" <<std::endl;
      std::cout<< "See LICENSE.txt for copyright info." <<std::endl <<std::endl;
      if (!parse.error() && argc!=0){
	valid_=1;
      }
    }

    bool is_valid() const {return valid_;}

    std::string filename() { //implicitly rather than explicitly const
      if (is_valid()){
	return std::string(parse.nonOption(0));
      }
      else {
	return std::string();
      }
    }

    long chi() const {
      if (is_valid() && options[CHI]){
	return stol(options[CHI].arg);
      }
      else {
	return 0;
      }
    }

    double trunc() const {
      if (is_valid() && options[TRUNC]){
	return stod(options[TRUNC].arg);
      }
      else {
	return 1.0e-16;
      }
    }

    void print(){
      if (!is_valid()) {
	option::printUsage(std::cout, usage_);
      }
      else {
	for (auto n=0; n<parse.nonOptionsCount();++n){
	  std::cout << "ARGUMENT " << n+1 << ": " << std::string(parse.nonOption(n)) <<std::endl;
	}
	for (auto&& o : options){
	  if (o)
	    for (option::Option* opt = o; opt; opt = opt->next()){
	      std::cout << "OPTION: " << std::string(opt->name,opt->namelen);
	      if (opt->arg)
		std::cout << " = " << opt->arg;
	      std::cout << std::endl;
	    }
	}
      }
    }
  };

  class iDMRG_Args : public Base_Args{
  private:
    unsigned long num_steps_;
    std::vector<short int> target_;

  public:
    iDMRG_Args(int argc, char* argv[]) : Base_Args(argc,argv,iDMRG_usage), num_steps_(0){
      
      //REQUIRE ONE NON OPTIONAL ARGUMENT, TREAT AS A FILENAME
      if (parse.nonOptionsCount()!=1 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
      else {
	valid_= (valid_==1);
      }
      if (is_valid()){
	if (options[NUMBER_OF_STEPS])
	  num_steps_=stoul(options[NUMBER_OF_STEPS].arg);
	if (options[TARGET]){
	  std::istringstream tss(options[TARGET].arg);
	  tss >> target_;
	}
      }
      print();
    }
    unsigned long number_of_steps() const {
      return num_steps_;
    }

    const std::vector<short int>& target() const {
      return target_;
    }

  };

  class fDMRG_Args : public Base_Args{
  private:
    unsigned int E_;
    unsigned int F_;
    double Weight_;
    unsigned int N_; //used by finite codes
    std::vector<short int> target_;

  public:
    fDMRG_Args(int argc, char* argv[]) : Base_Args(argc,argv,fDMRG_usage), E_(0), F_(0), Weight_(100.0),N_(0) {
      //REQUIRE TWO NON OPTIONAL ARGUMENTS, TREAT AS A FILENAME AND NUMBER OF VERTICES
      if (parse.nonOptionsCount()!=2 || std::string(parse.nonOption(0))==std::string("-") || options[NUMBER_OF_EXCITED].count()>1 ){
	std::cout << "Incorrect command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
      else {
	N_=stoul(parse.nonOption(1));
	valid_=(valid_==1);
      }
      //
      if (is_valid()){
	if (N_==0 || (N_ % 2)){
	  std::cout << "Illegal number of vertices requested: " << N_ << std::endl;
	  std::cout << "Must be a positive even value!" <<std::endl<<std::endl;
	  valid_=0;
	}
	if (options[WEIGHT_FACTOR]) {
	  E_=1; //at least 1 excited state to look for.
	  Weight_=stod(options[WEIGHT_FACTOR].arg);
	}
	if (options[NUMBER_OF_EXCITED])
	  E_=stoul(options[NUMBER_OF_EXCITED].arg);
	if (options[NUMBER_OF_SWEEPS])
	  F_=stoul(options[NUMBER_OF_SWEEPS].arg);
	if (options[TARGET]){
	  std::istringstream tss(options[TARGET].arg);
	  tss >> target_;
	}
      }
      print();
    };
    
    unsigned int num_vertices() const {return N_;}
    unsigned int num_excited() const {return E_;}
    unsigned int num_sweeps() const {return F_;}
    double weight_factor() const {return Weight_;}
    const std::vector<short int>& target() const {
      return target_;
    }

  };

  class iTEBD_Args : public Base_Args{
    unsigned long num_steps_;
    double step_size_;
    unsigned long trotter_order_;
    unsigned long measurement_interval_;
    std::string initial_unit_cell_;
    std::string c_number_filename_;

  public:
    iTEBD_Args(int argc, char* argv[]) : Base_Args(argc,argv,iTEBD_usage), num_steps_(1), step_size_(0.1), trotter_order_(2), measurement_interval_(1){
      
      //REQUIRE ONE NON OPTIONAL ARGUMENT, TREAT AS A FILENAME
      if (parse.nonOptionsCount()!=1 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
      else {
	valid_= (valid_==1);
      }
      if (is_valid()){
	if (options[NUMBER_OF_STEPS])
	  num_steps_=stoul(options[NUMBER_OF_STEPS].arg);
	if (options[TROTTER_ORDER])
	  trotter_order_=stoul(options[TROTTER_ORDER].arg);
	if (options[MEASUREMENT_INTERVAL])
	  measurement_interval_=stoul(options[MEASUREMENT_INTERVAL].arg);
	if (options[STEP_SIZE])
	  step_size_=stod(options[STEP_SIZE].arg);
	if (options[INITIAL_STATE_NAME])
	  initial_unit_cell_=std::string(options[INITIAL_STATE_NAME].arg);
	if (options[C_SPECIFIER]){
	  c_number_filename_=std::string(options[C_SPECIFIER].arg);
	}

	if (options[C_SPECIFIER] && options[INITIAL_STATE_NAME]){
	  valid_=0;
	  std::cout <<"Cannot specify initial state through options -c and -i simultaneously!"<<std::endl;
	}
      }
      print();
    }

    unsigned long number_of_steps() const {
      return num_steps_;
    }
    unsigned long trotter_order() const {
      return trotter_order_;
    }
    double step_size() const {
      return step_size_;
    }
    unsigned long measurement_interval() const {
      return measurement_interval_;
    }
    const std::string& initial_unit_cell() const {
      return initial_unit_cell_;
    }

    const std::string& c_number_filename() const {
      return c_number_filename_;
    }

  };

  class TEBD_Args : public Base_Args{
    unsigned long num_steps_;
    double step_size_;
    unsigned long trotter_order_;
    unsigned long measurement_interval_;
    std::string initial_state_name_;
    unsigned int N_; //used by finite codes
    std::vector<StringIndexPairs> finite_measurements_;
    std::string c_number_filename_;
    bool save_all_;

  public:
    TEBD_Args(int argc, char* argv[]) : Base_Args(argc,argv,TEBD_usage), num_steps_(1), step_size_(0.1), trotter_order_(2), measurement_interval_(1),N_(0),save_all_(0){
      
     if (parse.nonOptionsCount()!=2 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
      else {
	N_=stoul(parse.nonOption(1));
	valid_=(valid_==1);
      }
      //
      if (is_valid()){
	if (N_==0 || (N_ % 2)){
	  std::cout << "Illegal number of vertices requested: " << N_ << std::endl;
	  std::cout << "Must be a positive even value!" <<std::endl<<std::endl;
	  valid_=0;
	}
	if (options[NUMBER_OF_STEPS])
	  num_steps_=stoul(options[NUMBER_OF_STEPS].arg);
	if (options[TROTTER_ORDER])
	  trotter_order_=stoul(options[TROTTER_ORDER].arg);
	if (options[MEASUREMENT_INTERVAL])
	  measurement_interval_=stod(options[MEASUREMENT_INTERVAL].arg);
	if (options[STEP_SIZE])
	  step_size_=stod(options[STEP_SIZE].arg);
	if (options[INITIAL_STATE_NAME])
	  initial_state_name_=std::string(options[INITIAL_STATE_NAME].arg);
	if (options[FINITE_MEASUREMENT]){
	  for (option::Option* opt = options[FINITE_MEASUREMENT]; opt; opt = opt->next()){
	    StringIndexPairs temp;
	    std::istringstream ss(opt->arg);
	    ss >> temp;
	    //check all locations specified in temp
	    for (auto&& l : temp){
	      if (l.second < 1 || l.second >N_) {
		std::cout << "Specified measurement vertex " << l.second << " is outside bounds 1:" <<N_<<std::endl<<std::endl;
		valid_=0;
	      }
	    }
	    if (valid_)
	      finite_measurements_.emplace_back(temp);
	  }
	}
	
	if (options[C_SPECIFIER]){
	  c_number_filename_=std::string(options[C_SPECIFIER].arg);
	}

	if (options[C_SPECIFIER] && options[INITIAL_STATE_NAME]){
	  valid_=0;
	  std::cout <<"Cannot specify initial state through options -c and -i simultaneously!"<<std::endl;
	}
	
	if (options[SAVE_ALL]) {
	  save_all_=1; //set true
	  if (measurement_interval_!=1){
	    std::cout <<"Cannot specify both 'save-all' and a 'measurement-interval' not equal to 1." <<std::endl;
	    valid_=0;
	  }
	  else { 
	    measurement_interval_=1;
	  }
	}
      }
      print();
    }

    unsigned int num_vertices() const {return N_;}

    unsigned long number_of_steps() const {
      return num_steps_;
    }
    unsigned long trotter_order() const {
      return trotter_order_;
    }
    double step_size() const {
      return step_size_;
    }
    unsigned long measurement_interval() const {
      return measurement_interval_;
    }

    const std::string& initial_state_name() const {
      return initial_state_name_;
    }

    const std::vector<StringIndexPairs>& finite_measurements() const {
      return finite_measurements_;
    }

    const std::string& c_number_filename() const {
      return c_number_filename_;
    }

    bool save_all_flag() const {
      return save_all_;
    }

  };

  class iMEAS_Args : public Base_Args{
    unsigned long separation_;
    unsigned long nev_;
    bool use_filename_index_;
    bool two_point_;
    bool calc_entanglement_;
    unsigned long vert_entanglement_;
    std::vector<std::string> operator_filenames_;
    std::vector<std::string> files_;
    std::vector<short int> target_;
    std::string timefile_;
    std::string modelfile_;

  public:
    iMEAS_Args(int argc, char* argv[]) : Base_Args(argc,argv,iMEAS_usage),separation_(0),nev_(0),use_filename_index_(1),two_point_(0),calc_entanglement_(0),vert_entanglement_(0),timefile_(std::string()){
      if (valid_){
	if (parse.nonOptionsCount()<1 || std::string(parse.nonOption(0))==std::string("-")) {std::cout << "No files to process?" << std::endl<<std::endl; valid_=0;}
	else {
	  for (size_t f=0;f<parse.nonOptionsCount();++f){
	    files_.emplace_back(parse.nonOption(f));
	  }
	  if (options[SEPARATION]){ //has any separation been explicitly defined, even if zero? Then definitely a two point function
	    separation_=stoul(options[SEPARATION].arg);
	    two_point_=1;
	  }
	  if (options[NEV] && !options[TARGET])
	    valid_=0;
	  if (options[TARGET]){
	    std::istringstream tss(options[TARGET].arg);
	    tss >> target_;
	    nev_=1;
	    if (options[NEV])
	      nev_=stoul(options[NEV].arg);
	  }
	  if (options[NOINDEX])
	    use_filename_index_=0;
	  if (options[OPERATORFILE])
	    for (option::Option* opt = options[OPERATORFILE]; opt; opt = opt->next()){
	      operator_filenames_.emplace_back(opt->arg);
	    }
	  if (operator_filenames_.size()>1){//defined more than one op? then two point, on same vertex
	    two_point_=1;
	  }
	  if (options[ENTANGLEMENT])
	    calc_entanglement_=1;
	  if (options[VERTEX_ENTANGLEMENT])
	    vert_entanglement_=stoul(options[VERTEX_ENTANGLEMENT].arg);
	  if (options[TIMEFILE])
	    timefile_=std::string(options[TIMEFILE].arg);
	  if (options[IENERGY])
	    modelfile_=std::string(options[IENERGY].arg);
	}
      }
      else valid_=0;
      print();
    }

    unsigned long separation() const {return separation_;}
    const std::vector<std::string>& files() const {return files_;}
    const std::vector<std::string>& operator_filenames() const {return operator_filenames_;}
    bool use_filename_index() const {return use_filename_index_;}
    bool two_point() const {return two_point_;} //has a two point function been requested?
    unsigned long nev() const {return nev_;}
    bool calc_entanglement() const {return calc_entanglement_;}
    unsigned long vert_entanglement() const {return vert_entanglement_;}
    const std::vector<short int>& target() const {return target_;}
    const std::string& time_filename() const {return timefile_;}
    const std::string& model_filename() const {return modelfile_;}
  };

class fMEAS_Args : public Base_Args{
  unsigned int N_; //used by finite codes
  std::vector<StringIndexPairs> finite_measurements_;
  std::vector<std::string> state_names_;
  bool fdmrg_mode_;

  public:
    fMEAS_Args(int argc, char* argv[]) : Base_Args(argc,argv,fMEAS_usage) {
      
      if (parse.nonOptionsCount()<2 || (parse.nonOptionsCount()==2 && !options[FDMRG_MODE]) || (parse.nonOptionsCount()>=3 && options[FDMRG_MODE]) || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
     else {
       N_=stoul(parse.nonOption(1));
       if (!options[FDMRG_MODE]) {
	fdmrg_mode_=0;
	 for (size_t f=2;f<parse.nonOptionsCount();++f){
	   state_names_.emplace_back(parse.nonOption(f));
	 }
       }
      else
	fdmrg_mode_=1;
      
       valid_=(valid_==1);
     }
      //
     if (is_valid()){
       if (N_==0 || (N_ % 2)){
	 std::cout << "Illegal number of vertices requested: " << N_ << std::endl;
	 std::cout << "Must be a positive even value!" <<std::endl<<std::endl;
	 valid_=0;
       }
       if (options[FINITE_MEASUREMENT]){
	 for (option::Option* opt = options[FINITE_MEASUREMENT]; opt; opt = opt->next()){
	   StringIndexPairs temp;
	   std::istringstream ss(opt->arg);
	   ss >> temp;
	   //check all locations specified in temp
	   for (auto&& l : temp){
	     if (l.second < 1 || l.second >N_) {
	       std::cout << "Specified measurement vertex " << l.second << " is outside bounds 1:" <<N_<<std::endl<<std::endl;
	       valid_=0;
	     }
	   }
	   if (valid_)
	     finite_measurements_.emplace_back(temp);
	 }
       }
     }
     print();
    }

    unsigned int num_vertices() const {return N_;}
    const std::vector<std::string>& state_names() const {return state_names_;}
    const std::vector<StringIndexPairs>& finite_measurements() const {
      return finite_measurements_;
    }
  bool fdmrg_mode() const {return fdmrg_mode_;}

  };
  
  class Store_Args : public Base_Args{
  public:
    Store_Args(int argc, char* argv[]) : Base_Args(argc,argv,store_usage){
      if (parse.nonOptionsCount()!=1 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
      print();
    }
  };

  class TwoVE_Args : public Base_Args{
    unsigned long num_steps_;
    double step_size_;
    std::string initial_state_name_;
    std::vector<StringIndexPairs> finite_measurements_;
    std::string c_number_filename_;
    
  public:
    TwoVE_Args(int argc, char* argv[]) : Base_Args(argc,argv,TwoVE_usage), num_steps_(1), step_size_(0.1) {
      
      //Need at least the model specified on command line
     if (parse.nonOptionsCount()!=1 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
      else {
	valid_=(valid_==1);
      }
      //
      if (is_valid()){
	if (options[NUMBER_OF_STEPS])
	  num_steps_=stoul(options[NUMBER_OF_STEPS].arg);
	if (options[STEP_SIZE])
	  step_size_=stod(options[STEP_SIZE].arg);
	if (options[INITIAL_STATE_NAME])
	  initial_state_name_=std::string(options[INITIAL_STATE_NAME].arg);
	if (options[FINITE_MEASUREMENT]){
	  for (option::Option* opt = options[FINITE_MEASUREMENT]; opt; opt = opt->next()){
	    StringIndexPairs temp;
	    std::istringstream ss(opt->arg);
	    ss >> temp;
	    //check all locations specified in temp
	    for (auto&& l : temp){
	      if (l.second < 1 || l.second >2) {
		std::cout << "Specified measurement vertex " << l.second << " is not 1 or 2!" <<std::endl;
		valid_=0;
	      }
	    }
	    if (valid_)
	      finite_measurements_.emplace_back(temp);
	  }
	}
	if (options[C_SPECIFIER]){
	  c_number_filename_=std::string(options[C_SPECIFIER].arg);
	}

	if (options[C_SPECIFIER] && options[INITIAL_STATE_NAME]){
	  valid_=0;
	  std::cout <<"Cannot specify initial state through options -c and -i simultaneously!"<<std::endl;
	}
	  

      }
      print();
    }

    unsigned int num_vertices() const {return 2;}

    unsigned long number_of_steps() const {
      return num_steps_;
    }
    double step_size() const {
      return step_size_;
    }

    const std::string& initial_state_name() const {
      return initial_state_name_;
    }

    const std::vector<StringIndexPairs>& finite_measurements() const {
      return finite_measurements_;
    }

    const std::string& c_number_filename() const {
      return c_number_filename_;
    }

  };

}
#endif
