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
#include <utility>

#include "ajaj_common.hpp"

namespace ajaj {
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
    uMPXInt Step;
  public:
    const std::string filename;
    DataOutput(const std::string& f) : Step(0), filename(f) {outfile.open(filename.c_str(),ios::out | ios::trunc);outfile << setprecision(16);}
    ~DataOutput() {outfile.close();}
    bool push(const Data& results){
      if (outfile.is_open()){
	outfile << Step++;
	for (std::vector<double>::const_iterator cit=results.Real_measurements.begin();cit!=results.Real_measurements.end();++cit){
	  outfile << " " << *cit;
	}
	for (std::vector<std::complex<double> >::const_iterator cit=results.Complex_measurements.begin();cit!=results.Complex_measurements.end();++cit){
	  outfile << " " << real(*cit) << " " << imag(*cit) << " ";
	}
	outfile << std::endl;
	return 0;
      }
      else {
	return 1;
      }
    }
  };
}
#endif
