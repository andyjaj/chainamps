/** @file DMRG_routines.hpp
 * Driver functions for DMRG
 * Currently only two site DMRG is implemented.
 *
 */
#ifndef DMRG_H
#define DMRG_H

#include <vector>
#include <utility>
#include <string>
#include <sstream>


#include "ajaj_common.hpp"
#include "MPX.hpp"
#include "sparse_interface.hpp"
#include "data.hpp"

namespace ajaj {
  struct Prediction;
  class BlocksStructure;
  class SuperBlock;
  class ProjectorBlocks;
  //class DMRGBase;
  class iDMRG;
  class FiniteDMRG;
  class FiniteExcitedStates;
  class TwoVertexComponents;

  typedef std::pair<MPX_matrix,double> TensorWeightPair;

  /** Structure for state prediction vector. Used to speed up eigensolver.*/
  struct Prediction {
  public:
    SparseMatrix Guess; /**< The actual prediction vector */
    SparseMatrix LambdaL;
    SparseMatrix LambdaR; /**< Used when checking the overlap of the prediction vector (fidelity) with a calculated wavefunction.*/
  };

  /** A container of sorts, that stores blocks and provides interface to retrieve them.*/
  class BlocksStructure {
  private:
    std::string Name_; //name for file storage
    const EigenStateArray* SpectrumPtr_; //shouldn't change. Using a ptr to a const means we can use the default equals provided by the compiler
    MPX_matrix LeftBlock; //local storage
    MPX_matrix RightBlock;
    uMPXInt num_vertices_;
    uMPXInt numLeft_;
    uMPXInt numMiddle_;

  protected:
    bool save_left_block(uMPXInt l);
    bool save_right_block(uMPXInt r);
    bool save_left_block();
    bool save_right_block();
    bool save_blocks();
    void load_left_block();
    void load_right_block();
    std::string ResetName(const std::string& NewName, uMPXInt a){
      std::stringstream namestream;
      namestream << NewName << "_" << a;
      std::string OldName=Name_;
      Name_=namestream.str();
      return OldName;
    }

  public:
    BlocksStructure(const std::string&  Name, const EigenStateArray& Spectrum, uMPXInt num_vertices, uMPXInt numLeft, uMPXInt numMiddle) : Name_(Name), SpectrumPtr_(&Spectrum),LeftBlock(Spectrum),RightBlock(Spectrum),num_vertices_(num_vertices),numLeft_(numLeft),numMiddle_(numMiddle) {}
    BlocksStructure(const std::string&  Name, const EigenStateArray& Spectrum) : BlocksStructure(Name,Spectrum,0,0,0) {}

    uMPXInt size() const {return num_vertices_;}
    uMPXInt left_size() const {return numLeft_;}
    uMPXInt middle_size() const {return numMiddle_;}
    uMPXInt right_size() const {return num_vertices_-numLeft_-numMiddle_;}
    void ends(MPX_matrix&& leftend, MPX_matrix&& rightend);
    void initial_2(MPX_matrix&& newleft, MPX_matrix&& newright);
    void insert_2();
    void insert_2(MPX_matrix&& newleft, MPX_matrix&& newright);
    void set_2();
    virtual void shift_left(MPX_matrix&& newright);
    virtual void shift_right(MPX_matrix&& newleft);
    void shift_left(const std::string& OtherLeftName, MPX_matrix&& newright);
    void shift_right(MPX_matrix&& newleft, const std::string& OtherRightName);
    void set_blocks(uMPXInt ls, MPX_matrix&& newleft, MPX_matrix&& newright); //if blocks are generated by an algorithm other than grow

    const MPX_matrix& getLeftBlock() const {return LeftBlock;}
    const MPX_matrix& getRightBlock() const {return RightBlock;}
    const std::string& getName() const {return Name_;}
    const EigenStateArray& getSpectrum() const {return *SpectrumPtr_;}
  };

  /** Holds a superblock and can grow it or sweep it using DMRG methods. */
  class SuperBlock : public BlocksStructure{
  private:
    //For left and right blocks to be meaningful, there must be a Hamiltonian
    const MPO_matrix* H_ptr_;
    State TargetState_;
    std::string DensityFileName_;

  protected:
    std::vector<double> previous_lambda_;
    double fidelity_;
    Prediction pred_;
  public:
    MPSDecomposition CentralDecomposition;

    SuperBlock(const std::string& Name, const MPO_matrix& H, const State& TargetState) : BlocksStructure(Name,H.GetPhysicalSpectrum(),0,0,0),H_ptr_(&H),TargetState_(TargetState),CentralDecomposition(H.getPhysicalSpectrum()) {
      std::stringstream dnamestream;
      dnamestream << getName() << "_One_Vertex_Densities.dat";
      DensityFileName_=dnamestream.str();
      std::ofstream DensityFileStream_;
      DensityFileStream_.open(DensityFileName_.c_str(),ios::out | ios::trunc);
      for (auto&& i : H.getPhysicalSpectrum().Energies()){
	DensityFileStream_ << i << " ";
      }
      DensityFileStream_ << std::endl;
      DensityFileStream_.close();
    }
    const MPO_matrix& getH() const {return *H_ptr_;}
    const MPSDecomposition& getCentralDecomposition() const {return CentralDecomposition;}
    const std::vector<double>& getPreviousLambda() const {return previous_lambda_;}
    const State& getTargetState() const {return TargetState_;}
    void push_density() const;

    double getTruncation() const {return CentralDecomposition.Truncation;}

    virtual Data initialise(uMPXInt chi=0, double smin=0.0); //two vertex initialisation, returns some two vertex measurements
    virtual Data grow_two_vertex(uMPXInt chi=0, double smin=0.0);
    virtual Data move_right_two_vertex(uMPXInt chi=0, double smin=0.0);
    virtual Data move_left_two_vertex(uMPXInt chi=0, double smin=0.0);
  };

  class ProjectorBlocks : public BlocksStructure{
  private:
    double Weight_;
    const std::string ProjectorStateName_;
    std::pair<MPS_matrix,MPS_matrix> ProjectorStatePair_;
    uMPXInt left_mark_;
    uMPXInt right_mark_;

    std::pair<MPS_matrix,MPS_matrix> FetchProjectorStatePair(uMPXInt ls);
    std::string MakeName(const std::string& Name) const{
      std::stringstream namestream;
      namestream << "Projector_" << Name;
      return namestream.str();
    }
  public:
    ProjectorBlocks(const std::string&  Name, const EigenStateArray& Spectrum, uMPXInt num_vertices, uMPXInt numLeft, uMPXInt numMiddle, double Weight) : BlocksStructure(MakeName(Name),Spectrum,num_vertices,numLeft,numMiddle),Weight_(Weight),ProjectorStateName_(Name),ProjectorStatePair_(FetchProjectorStatePair(numLeft)),left_mark_(numLeft),right_mark_(num_vertices-numMiddle-numLeft) {
      //projector state pair knows to incorporate the singular values at the centre of the system
      const MPS_matrix& psi_f=ProjectorStatePair_.first;
      const MPS_matrix& psi_s=ProjectorStatePair_.second;
      set_blocks(left_size(),MPX_matrix(Spectrum,psi_f.getInwardMatrixIndex(),std::vector<double>(psi_f.getInwardMatrixIndex().size(),1.0)),std::move(MPX_matrix(Spectrum,psi_s.getOutwardMatrixIndex(),std::vector<double>(psi_s.getOutwardMatrixIndex().size(),1.0)).Transpose()));
    }
    double getWeight() const {return Weight_;}
    TensorWeightPair makePTensor() const;
    void move_left_two_vertex(const MPS_matrix& RightMatrix);
    void move_right_two_vertex(const MPS_matrix& LeftMatrix);
  };

  class iDMRG : public SuperBlock {
  private:
    DataOutput& output_ref_;
  public:
    iDMRG(const std::string& Name, const MPO_matrix& H, const State& TargetState, DataOutput& resultsref) : SuperBlock(Name,H,TargetState),output_ref_(resultsref) {};
    void run(uMPXInt number_of_steps=0, double convergence_criterion=0.0,  uMPXInt chi=0, double smin=0.0); /**< Perform infinite algorithm growth steps*/
    double fidelity() const {return fidelity_;}
  };

  class FiniteDMRG : public SuperBlock {
  private:
    double chi_;
    double smin_;
    DataOutput& output_ref_;
  public:
    FiniteDMRG(SuperBlock& previous, DataOutput& resultsref) : SuperBlock(std::move(previous)),output_ref_(resultsref) {}
    FiniteDMRG(iDMRG& infrun, DataOutput& resultsref) : SuperBlock(std::move(infrun)),output_ref_(resultsref) {}
    void run(uMPXInt num_sweeps, uMPXInt chi=0, double smin=0.0);
    double chi() const {return chi_;}
    double smin() const {return smin_;}
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
    ExcitedStateFiniteDMRG(const std::string& Name, FiniteDMRG& FD, double Weight, DataOutput& resultsref) : SuperBlock(std::move(FD)),output_ref_(resultsref),BaseName_(Name),GSName_(ResetName(BaseName_,1)),HBlocksName_(GSName_),init_flag_(1),PBlocks_(1,ProjectorBlocks(HBlocksName_,getSpectrum(),size(),left_size(),middle_size(),Weight)) {}
    void init(double chi, double smin);
    void run(uMPXInt number_of_sweeps, uMPXInt chi, double smin);
    Data move_right_two_vertex(uMPXInt chi=0, double smin=0.0);
    Data move_left_two_vertex(uMPXInt chi=0, double smin=0.0);

    void next_state(double Weight){
      //push_back projectorblocks       //update name
      HBlocksName_=getName();
      PBlocks_.push_back(ProjectorBlocks(ResetName(BaseName_,PBlocks_.size()+1),getSpectrum(),size(),left_size(),middle_size(),Weight));
      init_flag_=1;
    };

  };

  /** Form a left end, open boundary conditions, Hamiltonian MPO_matrix. The leftmost matrix index will be a dummy. */
  inline MPO_matrix LeftOpenBCHamiltonian(const MPO_matrix& H){return H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(1,H.dimsvector()[1]-1)));}
  /** Form a right end, open boundary conditions, Hamiltonian MPO_matrix. The rightmost matrix index will be a dummy. */
  inline MPO_matrix RightOpenBCHamiltonian(const MPO_matrix& H){return H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(3,0)));}
  /** Do an SVD, with optional truncation, on an MPX_matrix object (with two physical indices and two matrix indices), forming two new MPS_matrix objects.*/
  inline MPSDecomposition TwoVertexSVD(const MPX_matrix& M, size_t bond_dimension=0, double min_s_val=0.0){return MPSDecomposition(M.SVD(bond_dimension,min_s_val));}
  /** Special initial left block with dummy indices, used for first DMRG step.*/
  MPX_matrix MakeInitialLeftBlock(const MPO_matrix& LeftH, const MPS_matrix& A);
  MPX_matrix MakeInitialRightBlock(const MPO_matrix& RightH, const MPS_matrix& B);
  /** Forms the two vertex Hamiltonian */
  MPX_matrix TwoVertexInitialHamiltonian(const MPO_matrix& LeftH, const MPO_matrix& RightH);
  /** Solves the two vertex Hamiltonian eigenproblem and makes a two vertex wavefunction, output requires SVD. */
  MPX_matrix TwoVertexInitialWavefunction(const MPO_matrix& LeftH, const MPO_matrix& RightH, const State& TargetSector, Data& result);
  /** Makes a prediction vector for the eigensolver*/
  Prediction MakePrediction(const MPSDecomposition& Decomp, const std::vector<double>& PreviousLambda);
  /** Updates the left block*/
  MPX_matrix MakeLeftBlock(const MPX_matrix& LB, const MPX_matrix& H, const MPS_matrix& A);
  MPX_matrix MakeRightBlock(const MPX_matrix& RB, const MPX_matrix& H, const MPS_matrix& B);
  MPX_matrix MakeDummyLeftBlock(const MPO_matrix& H, const State& TargetState);
  MPX_matrix MakeDummyRightBlock(const MPO_matrix& H, const State& TargetState);
  /** Solves the superblock Hamiltonian eigen problem, output requires SVD. */
  MPX_matrix TwoVertexWavefunction(const MPX_matrix& LeftBlock, const MPO_matrix& H, const MPX_matrix& RightBlock, const std::vector<ProjectorBlocks>* ProjectorBlocksPtr, MPXInt NumVertices, Data& result, SparseMatrix* guessptr=nullptr);
  /** checks overlap between prediction vector and new wavefunction */
  double CheckConvergence(const Prediction& guess,const std::vector<double>& Lambda);

  class TwoVertexComponents{
  public:
    const MPX_matrix& LeftBlock;
    const MPO_matrix& H;
    const MPX_matrix& RightBlock;
    const std::vector<ProjectorBlocks>* ProjectorsPtr;//needs to numbers for left and right in order to load
    const State* TargetStatePtr;

    std::vector<MPXIndex> indices;
    std::vector<MPXInt> allowed_indices;
    std::vector<MPXInt> left_allowed_indices;
    std::vector<std::array<MPXInt,2> > rows_and_cols;
    uMPXInt vrows;
    uMPXInt vcols;
    MPX_matrix LeftPart;
    MPX_matrix RightPart;
    std::vector<TensorWeightPair> ProjectorTensors;

    TwoVertexComponents(const MPX_matrix& L, const MPO_matrix& HMPO, const MPX_matrix& R, const std::vector<ProjectorBlocks>* P=nullptr, const State* StatePtr=nullptr);
    //TwoVertexComponents(const MPX_matrix& L, const MPO_matrix& HMPO, const MPX_matrix& R, const State* StatePtr) : TwoVertexComponents(L,HMPO,R,nullptr,StatePtr) {}


    SparseHED HED(MPXInt numevals, char which[3],const SparseMatrix* initial=NULL) const;
    MPXInt length()const {return m_length;}
  private:
    uMPXInt m_length;
  };

  void TwoVertexMPOMPSMultiply(const TwoVertexComponents* array, std::complex<double> *in, std::complex<double> *out);

}

#endif
