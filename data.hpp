/** @file data.hpp
 * structure for outputting data to file
 *
 */

#ifndef DATA_H
#define DATA_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>
#include <sstream>
#include <utility>

#include "ajaj_common.hpp"

namespace ajaj {

  inline std::string StripName(const std::string& fullpath){
    return fullpath.find("/")!=std::string::npos ?  fullpath.substr(fullpath.rfind("/")+1,fullpath.length()) : fullpath;
  }

  inline std::string OutputName(const std::string& InfileName, const std::string& ResultTypeName){
    std::ostringstream outss;
    outss << StripName(InfileName) <<"_"<<ResultTypeName;
    return outss.str();
  }

  struct Data{
    std::vector<double> Real_measurements;
    std::vector<std::complex<double> > Complex_measurements;
    Data(){}
    Data(const std::vector<double>& r) {Real_measurements=r;Complex_measurements.resize(0);}
    Data(const std::vector<std::complex<double> >& c)  {Real_measurements.resize(0);Complex_measurements=c;}
    Data(const std::vector<double>& r, const std::vector<std::complex<double> >& c) {Real_measurements=r;Complex_measurements=c;}
    Data(double r) {Real_measurements=std::vector<double>(1,r);Complex_measurements.resize(0);}
    Data(std::complex<double> c)  {Real_measurements.resize(0);Complex_measurements=std::vector<std::complex<double> >(1,c);}
    Data(double r, std::complex<double> c)   {Real_measurements=std::vector<double>(1,r);Complex_measurements=std::vector<std::complex<double> >(1,c);}
    void swap(Data& other){
      std::swap(Real_measurements,other.Real_measurements);
      std::swap(Complex_measurements,other.Complex_measurements);
    }
    Data& operator=(Data other){
      swap(other);
      return *this;
    }
  };

  class DataOutput{
  private:
    std::ofstream outfile;
    uMPXInt Step; //for internal monitoring of the step
  public:
    const std::string filename;
    std::string comment;
    DataOutput(const std::string& f,const std::string& cs=std::string()) : Step(0), filename(f),comment(cs) {
      outfile.open(filename.c_str(),ios::out | ios::trunc);outfile << setprecision(16);
      if (!cs.empty()){
	outfile << "# "<< comment << std::endl;
      }
    }
    ~DataOutput() {outfile.close();}
    bool push(uMPXInt i,const Data& results){
      if (outfile.is_open()){
	outfile << i;
	for (std::vector<double>::const_iterator cit=results.Real_measurements.begin();cit!=results.Real_measurements.end();++cit){
	  outfile << " " << *cit;
	}
	for (std::vector<std::complex<double> >::const_iterator cit=results.Complex_measurements.begin();cit!=results.Complex_measurements.end();++cit){
	  outfile << " " << real(*cit) << " " << imag(*cit) << " ";
	}
	outfile << std::endl;
	Step++;
	return 0;
      }
      else {
	return 1;
      }
    }
    bool push(const Data& results){
      return push(Step,results);
    }
  };
}
#endif
