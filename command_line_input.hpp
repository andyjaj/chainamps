/** @file command_line_input.hpp
 * Initial input from command line can be given in the form of a spectrum and matrix elements
 * A model still needs to be defined though for the Hamiltonian to be formed correctly.
 */
#ifndef COMMAND_LINE_INPUT_H
#define COMMAND_LINE_INPUT_H
#include <iostream>
#include <string>
#include <vector>
#include "optionparser/optionparser.h" //The Lean Mean C++ Option Parser. See .h file for license info.

namespace ajaj {

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
    static option::ArgStatus PositiveDouble(const option::Option& option, bool msg)
    {
      char* endptr = nullptr;
      if (option.arg != 0 && stod(option.arg)>0){};
      if (endptr != option.arg && *endptr == 0)
	return option::ARG_OK;

      if (msg) std::cout << "Option '" << std::string(option.name,option.namelen) << "' requires a positive numeric argument" <<std::endl;
      return option::ARG_ILLEGAL;
    }
  };

  enum optionIndex {UNKNOWN,CHI,NUMBER_OF_STEPS,MINS,NUMBER_OF_EXCITED,NUMBER_OF_SWEEPS,WEIGHT_FACTOR,TROTTER_ORDER,TIME_STEPS,STEP_SIZE,MEASUREMENT_INTERVAL,INITIAL_STATE_NAME,SEPARATION,NOINDEX,OPERATORFILE};

  const option::Descriptor iDMRG_usage[4] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: iDMRG_DRV.bin [-B <number> -N <number>] <model_filename>"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, \t--bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
      {NUMBER_OF_STEPS,0,"N","number-of-steps",Arg::PositiveNumeric,"  -N <number>, \t--number-of-steps=<number>"
       "  \tThe number of infinite volume steps, >= 0" },
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor fDMRG_usage[6] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: fDMRG_DRV.bin [-B <number> -X <number> -F <number> -W <number>] <model_filename> <number of vertices/chains> \n  <number of vertices/chains> must be EVEN.\n"
	"Options:"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, \t--bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
      {NUMBER_OF_EXCITED,0,"X","excited-states",Arg::PositiveNumeric,"  -X <number>, \t--excited-states=<number>"
       "  \tThe number of excited states to calculate, >= 0. If this is specified but a projective weight factor is not, then a default value of 100.0 is used." },
      {NUMBER_OF_SWEEPS,0,"F","finite-size-sweeps",Arg::PositiveNumeric,"  -F <number>, \t--finite-size-sweeps=<number>"
       "  \tThe number of finite size sweeps to perform, >= 0" },
      {WEIGHT_FACTOR,0,"W","weight-factor",Arg::PositiveDouble,"  -W <number>, \t--weight-factor=<number>"
       "  \tThe weight factor if calculating excited states. Must be > 0. Specifying this indicates that the number of requested excited states is at least 1." },
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor iTEBD_usage[7] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: iTEBD_DRV.bin [-B <number> -n <number> -s <number> -O <number>] <model_filename>"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, \t--bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
      {NUMBER_OF_STEPS,0,"n","time-steps",Arg::PositiveNumeric,"  -n <number>, \t--time-steps=<number>"
       "  \tThe number of time steps. Default is 1." },
      {STEP_SIZE,0,"s","step-size",Arg::PositiveDouble,"  -s <number>, \t--step-size=<number>"
       "  \tThe step size. Default is 0.1" },
      {TROTTER_ORDER,0,"O","trotter-order",Arg::PositiveNumeric,"  -O <number>, \t--trotter-order=<number>"
       "  \tThe Trotter order (currently 1 or 2). Second order (2) is default." },
      {MEASUREMENT_INTERVAL,0,"m","measurement-interval",Arg::PositiveNumeric,"  -m <number>, \t--measurement-interval=<number>"
       "  \tMeasurement at every <number> steps. Default is 1 (measurement at every step)."},
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor TEBD_usage[8] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: TEBD_DRV.bin [-B <number> -n <number> -s <number> -O <number> -i <initial_state_name>] <model_filename> <number of vertices/chains> \n  <number of vertices/chains> must be EVEN.\n"},
      {CHI,0,"B","bond-dimension",Arg::PositiveNumeric,"  -B <number>, \t--bond-dimension=<number>"
       "  \tThe maximum bond dimension, >= 0. If 0, then ignored." },
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
      { 0, 0, 0, 0, 0, 0 }
    };

  const option::Descriptor iMEAS_usage[5] =
    {
      {UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: UNITCELL_MEASURE.bin [-D <number>] <unitcell_filename1> ... \n"},
      {OPERATORFILE,0,"O","operator filename",Arg::NonEmpty,"  -O <filename>, \t--operator-file=<filename>"
       "  \tDistance between two vertex measurements." },
      {SEPARATION,0,"S","separation",Arg::PositiveNumeric,"  -S <number>, \t--separation=<number>"
       "  \tDistance between two vertex measurements." },
      {NOINDEX,0,"X","No index",Arg::None,"  -X, \t--no-index"
       "  \tDon't extract index from filenames." },
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
      std::cout<< "See LICENSE.txt for copyright info." <<std::endl;
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
  public:
    unsigned long num_steps_;

    iDMRG_Args(int argc, char* argv[]) : Base_Args(argc,argv,iDMRG_usage), num_steps_(0){
      
      //REQUIRE ONE NON OPTIONAL ARGUMENT, TREAT AS A FILENAME
      if (parse.nonOptionsCount()!=1 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect number of command line arguments." << std::endl <<std::endl;
	valid_=0;
      }
      else {
	valid_= (valid_==1);
      }
      if (is_valid()){
	if (options[NUMBER_OF_STEPS])
	  num_steps_=stoul(options[NUMBER_OF_STEPS].arg);
      }
      print();
    }
    unsigned long number_of_steps() const {
      return num_steps_;
    }

  };

  class fDMRG_Args : public Base_Args{
  private:
    unsigned int E_;
    unsigned int F_;
    double Weight_;
    unsigned int N_; //used by finite codes


  public:
    fDMRG_Args(int argc, char* argv[]) : Base_Args(argc,argv,fDMRG_usage), E_(0), F_(0), Weight_(100.0),N_(0) {
      //REQUIRE TWO NON OPTIONAL ARGUMENTS, TREAT AS A FILENAME AND NUMBER OF VERTICES
      if (parse.nonOptionsCount()!=2 || std::string(parse.nonOption(0))==std::string("-") || options[NUMBER_OF_EXCITED].count()>1 ){
	std::cout << "Incorrect number of command line arguments." << std::endl <<std::endl;
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
	  std::cout << "Must be a positive even value!" <<std::endl;
	  valid_=0;
	}
	if (options[WEIGHT_FACTOR]) {
	  E_=1; //at least 1.
	  Weight_=stod(options[WEIGHT_FACTOR].arg);
	}
	if (options[NUMBER_OF_EXCITED])
	  E_=stoul(options[NUMBER_OF_EXCITED].arg);
	if (options[NUMBER_OF_SWEEPS])
	  F_=stoul(options[NUMBER_OF_SWEEPS].arg);
      }
      print();
    };
    
    unsigned int num_vertices() const {return N_;}
    unsigned int num_excited() const {return E_;}
    unsigned int num_sweeps() const {return F_;}
    double weight_factor() const {return Weight_;}

  };

  class iTEBD_Args : public Base_Args{
    unsigned long num_steps_;
    double step_size_;
    unsigned long trotter_order_;
    unsigned long measurement_interval_;

  public:
    iTEBD_Args(int argc, char* argv[]) : Base_Args(argc,argv,iTEBD_usage), num_steps_(1), step_size_(0.1), trotter_order_(2), measurement_interval_(1){
      
      //REQUIRE ONE NON OPTIONAL ARGUMENT, TREAT AS A FILENAME
      if (parse.nonOptionsCount()!=1 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect number of command line arguments." << std::endl <<std::endl;
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
	  measurement_interval_=stod(options[MEASUREMENT_INTERVAL].arg);
	if (options[STEP_SIZE])
	  step_size_=stod(options[STEP_SIZE].arg);
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
  };

  class TEBD_Args : public Base_Args{
    unsigned long num_steps_;
    double step_size_;
    unsigned long trotter_order_;
    unsigned long measurement_interval_;
    std::string initial_state_name_;
    unsigned int N_; //used by finite codes

  public:
    TEBD_Args(int argc, char* argv[]) : Base_Args(argc,argv,TEBD_usage), num_steps_(1), step_size_(0.1), trotter_order_(2), measurement_interval_(1),N_(0){
      
     if (parse.nonOptionsCount()!=2 || std::string(parse.nonOption(0))==std::string("-")){
	std::cout << "Incorrect number of command line arguments." << std::endl <<std::endl;
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
	  std::cout << "Must be a positive even value!" <<std::endl;
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

  };

  class iMEAS_Args : public Base_Args{
    unsigned long separation_;
    bool use_filename_index_;

    std::vector<std::string> operator_filenames_;
    std::vector<std::string> files_;

  public:
    iMEAS_Args(int argc, char* argv[]) : Base_Args(argc,argv,iMEAS_usage), separation_(0),use_filename_index_(1){
      if (is_valid()){
	for (size_t f=0;f<parse.nonOptionsCount();++f){
	  files_.emplace_back(parse.nonOption(f));
	}
	if (options[SEPARATION])
	  separation_=stoul(options[SEPARATION].arg);
	if (options[NOINDEX])
	  use_filename_index_=0;
	if (options[OPERATORFILE])
	  for (option::Option* opt = options[OPERATORFILE]; opt; opt = opt->next()){
	    operator_filenames_.emplace_back(opt->arg);
	  }
      }
      print();
    }

    unsigned long separation() const {return separation_;}
    const std::vector<std::string>& files() const {return files_;}
    const std::vector<std::string>& operator_filenames() const {return operator_filenames_;}
    bool use_filename_index() const {return use_filename_index_;}

  };

}
#endif
