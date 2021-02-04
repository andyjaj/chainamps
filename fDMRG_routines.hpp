#ifndef FDMRGROUTINES_H
#define FDMRGROUTINES_H

#include "DMRG_routines.hpp" //MPX_matrix etc.
#include "iDMRG_routines.hpp" //MPX_matrix etc.

namespace ajaj {
  
  class FiniteDMRG : public SuperBlock {
  private:
    double chi_;
    double truncation_;
    DataOutput& output_ref_;
  public:
    FiniteDMRG(SuperBlock& previous, DataOutput& resultsref) : SuperBlock(std::move(previous)),output_ref_(resultsref) {}
    FiniteDMRG(iDMRG& infrun, DataOutput& resultsref) : SuperBlock(std::move(infrun)),output_ref_(resultsref) {}
    
    void run(uMPXInt num_sweeps, uMPXInt chi=0, double truncation=0.0);
    double chi() const {return chi_;}
    double truncation() const {return truncation_;}
  };

  class ExcitedStateFiniteDMRG : public SuperBlock {
  private:
    DataOutput& output_ref_;
    std::string BaseName_;
    const std::string GSName_;
    std::string HBlocksName_; //need this to load H blocks
    bool init_flag_;
    std::vector<ProjectorBlocks> PBlocks_;
  public:
    //steal resources from finite dmrg object if possible
    ExcitedStateFiniteDMRG(const std::string& Name, FiniteDMRG& FD, double Weight, DataOutput& resultsref) : SuperBlock(std::move(FD)),output_ref_(resultsref),BaseName_(Name),GSName_(SetName(BaseName_,1)),HBlocksName_(GSName_),init_flag_(1),PBlocks_(1,ProjectorBlocks(HBlocksName_,getSpectrum(),size(),left_size(),middle_size(),Weight)) {}
    void init(double chi, double truncation, bool converge=1);
    void run(uMPXInt number_of_sweeps, uMPXInt chi, double truncation);
    Data move_right_two_vertex(uMPXInt chi=0, double truncation=0.0, bool converge=1);
    Data move_left_two_vertex(uMPXInt chi=0, double truncation=0.0, bool converge=1);

    void next_state(double Weight){
      //push_back projectorblocks       //update name
      HBlocksName_=getName();
      PBlocks_.push_back(ProjectorBlocks(SetName(BaseName_,PBlocks_.size()+1),getSpectrum(),size(),left_size(),middle_size(),Weight));
      init_flag_=1;
    };

  };

}
#endif
