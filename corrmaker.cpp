#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#include <string>
#include <regex>
#include <vector>
#include <array>
#include <algorithm>

#include "optionparser/optionparser.h" //The Lean Mean C++ Option Parser. See .h file for license info.

static const std::vector<std::pair<double,double> > momenta={{0.0,0.0},{M_PI,0.0},{0.0,M_PI},{M_PI,M_PI}};
static const size_t nummom(4);

static const size_t chain_length(4);

enum  optionIndex { UNKNOWN };
const option::Descriptor usage[] =
  {
    {UNKNOWN, 0,"" , ""    ,option::Arg::None, "USAGE: example [options]\n\n"},
    {0,0,0,0,0,0}
  };

struct datapoint {
  size_t Label;
  std::vector<std::complex<double> > Values;

  datapoint() {}

  datapoint(size_t l, const std::vector<std::complex<double> >& v) : Label(l),Values(v) {}
  datapoint(size_t l, const std::complex<double>& c) : Label(l),Values(std::vector<std::complex<double> >(1,c)){}
};

std::complex<double> phase(const std::pair<double,double>& Q, const std::pair<int,int>& coord){
  double factor(Q.first*coord.first+Q.second*coord.second);
  return exp(std::complex<double>(0.0,factor));
}

int main(int argc, char** argv){

  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats  stats(usage, argc, argv);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);

  option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);
  
  if (parse.error())
    return 1;

  //std::vector<ajaj::Data> results;
  std::vector<datapoint> results;
  std::vector<std::string> filenames;

  for (size_t f=0;f<parse.nonOptionsCount();++f){
    filenames.emplace_back(parse.nonOption(f));
  }

  std::regex x_pattern("_([0-9]+)\\.dat$"); //regex to match number before filename only.
  std::regex y_pattern("@1@([0-9]+)"); //regex to match number before filename only.

  bool first(1);
  size_t max_points(0);

  for (auto&& fn : filenames){
    //get an index

    unsigned long x(0);
    unsigned long y(0);

    //strip length data from file to get x,y pair

    std::smatch x_sm;
    if (std::regex_search(fn,x_sm,x_pattern) && x_sm.size()>1){
      x=std::stoul(x_sm.str(1).c_str());
    }

    std::smatch y_sm;
    if (std::regex_search(fn,y_sm,y_pattern) && y_sm.size()>1){
      y=std::stoul(y_sm.str(1).c_str());
    }

    double dblfactor=(y!=0 && y<chain_length/2) ? 2.0 : 1.0;

    std::cout << fn << std::endl;
    
    std::cout << x << " " << y <<std::endl;
    std::pair<int,int> coord(x,y);
    for (auto&& m: momenta){
      std::cout << phase(m,coord) << " ";      
    }
    std::cout <<std::endl;
    if (y!=0 && y<chain_length/2){
      std::cout << "Needs doubling" <<std::endl;
    }

    //open file
    std::ifstream infile;
    infile.open(fn.c_str(),std::ios::in);
    //parse file and save
    if (infile.is_open()){
      std::string comment;
      getline(infile,comment);
      //skip first line of file
      std::vector<datapoint> filedata;
      unsigned long index;
      double repart;
      double impart;
      std::string dataline;
      unsigned long num_lines(0);
      while (getline(infile,dataline)){
	++num_lines;
	std::istringstream iss(dataline);
	if (iss >> index >> repart >>impart){
	  filedata.emplace_back(index,std::complex<double>(repart,impart));
	}
	else {
	  std::cout << "Data corrupt!" <<std::endl;
	  exit(1);
	}
      }
      std::cout << filedata.size() <<std::endl;

      if (first){
	max_points=filedata.size();
	for (auto&& fd : filedata){
	  results.emplace_back(datapoint(fd.Label,std::vector<std::complex<double> >()));
	  //emplace_back results
	  for (auto&& m : momenta){
	    results.back().Values.push_back(dblfactor*phase(m,coord)*fd.Values[0]);
	  }
	}
	first=0;
      }
      else {
	if (filedata.size()<max_points) {max_points=filedata.size(); results.resize(max_points);}
	for (size_t i=0;i<max_points;++i){
	  for (size_t j=0;j<momenta.size();++j){
	    results[i].Values[j]+=dblfactor*phase(momenta[j],coord)*filedata[i].Values[0];
	    
	  }
	}
      }   
    }
  }

  //write out

  std::ofstream outfile;
  outfile.open("corr.dat",std::ios::out | std::ios::trunc);
  if (outfile.is_open()){
    for (auto&& r : results){
      outfile << r.Label << " ";
      for (auto&& v :r.Values){
	outfile << v.real() <<" " << v.imag() << " ";
      }
      outfile << std::endl;
    }
  }
  //store results

}
