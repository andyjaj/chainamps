#ifndef OLD_XXX_H
#define OLD_XXX_H

#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>
#include <sstream>
#include <limits>

//#include "../common_defs.hpp"
//#include "../sparse_interface.hpp"
//#include "../vertex.hpp"
//#include "../MPX.hpp"

namespace oldxxx {
  ajaj::MPO_matrix MakeHamiltonian(const ajaj::Vertex& modelvertex, const ajaj::CouplingArray& couplingparams){
    //Lower triangular MPO
    // I        0  0  0  0
    // JxyS+/2  0  0  0  0
    // JxyS-/2  0  0  0  0
    // JzSz     0  0  0  0
    // HV       S- S+ Sz I        //note the charges Q, are such that Q[Psidagger[i]]+Q[Psi'[i]]=0

    ajaj::QNCombinations differencecombinations(modelvertex.Spectrum,1); //1 means use difference

    ajaj::MPXInt lineardim(modelvertex.Spectrum.size()*(2+3*differencecombinations.size()));
    ajaj::MPXInt offset_to_last_block=modelvertex.Spectrum.size()*(1+3*differencecombinations.size()); //offset to get to the last row of operators
    ajaj::SparseMatrix M(lineardim,lineardim,lineardim);

    //start with the really easy bits, the Identities, I, and the vertex Hamiltonian HV
    for (size_t i=0;i<modelvertex.Spectrum.size();++i){
      M.entry(i,i,1.0);
      M.entry(i+offset_to_last_block,i,modelvertex.Spectrum[i].en);
      M.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
    }

    //now do the non trivial parts
    double Nxxx=modelvertex.Parameters[0].Value;
    //J1=Jxy=Jz
    
    for (auto&& c : couplingparams){
      if (c.Value.imag()!=0.0){
	std::cout << "ERROR: coupling param " << c << " must have a real value. Imag part=" << c.Value.imag() <<std::endl;
      }
    }
    double Jz=couplingparams[0].Value.real()*Nxxx;
    double halfJxy=couplingparams.size()>1 ? couplingparams[1].Value.real()*0.5*Nxxx : Jz*0.5;
    double J2z=couplingparams.size()>2 ? couplingparams[2].Value.real()*Nxxx : 0.0;
    double halfJ2xy=couplingparams.size()>3 ? couplingparams[3].Value.real()*0.5*Nxxx : J2z*0.5;
    std::cout << Jz/Nxxx << " " << 2.0*halfJxy/Nxxx << " " << J2z/Nxxx << " " << 2.0*halfJ2xy/Nxxx << std::endl;

    //Just use S+ and its conjugate?
    ajaj::Sparseint Sm_col_offset=1; //+1 for identity matrix in first block
    ajaj::Sparseint Sm_row_offset=differencecombinations.size()+1;
    ajaj::Sparseint Sp_col_offset=differencecombinations.size()+1; //+1 for identity matrix in first block
    ajaj::Sparseint Sp_row_offset=1; //reversed order
    for (ajaj::Sparseint col=0;col<modelvertex.Operators[0].MatrixElements.cols();++col){
      for (ajaj::Sparseint p=modelvertex.Operators[0].MatrixElements.get_p(col);p<modelvertex.Operators[0].MatrixElements.get_p(col+1);++p){
	ajaj::Sparseint i=modelvertex.Operators[0].MatrixElements.get_i(p);
	ajaj::State diffstate=modelvertex.Spectrum[i]-modelvertex.Spectrum[col];
	//now look through charge groups to find MPO indices
	ajaj::Sparseint MPO_subcol=0;
	ajaj::Sparseint MPO_subrow=0;
	ajaj::Sparseint MPO_subcolSp=0;
	ajaj::Sparseint MPO_subrowSp=0;
	if (diffstate==-diffstate){ //involution block
	  for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	    if (diffstate==differencecombinations.InvolutionPairs[l].PairState){
	      MPO_subcol=l;
	      MPO_subrow=l; //involution pairs don't get reversed
	      MPO_subcolSp=l;
	      MPO_subrowSp=l;
	      break;
	    }
	  }
	}
	else {
	  for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	    if (diffstate==differencecombinations.OrderedPairs[l].PairState){
	      MPO_subcol=differencecombinations.InvolutionPairs.size()+l;
	      MPO_subrow=differencecombinations.size()-l-1;
	      MPO_subcolSp=MPO_subrow;
	      MPO_subrowSp=MPO_subcol;
	      break;
	    }
	  }
	}
	std::complex<double> x=modelvertex.Operators[0].MatrixElements.get_x(p);
	double momentumfactor=2.0*halfJ2xy*cos(2.0*M_PI*diffstate[1]/Nxxx);
	//lowest block row
	M.entry(offset_to_last_block+i,modelvertex.Spectrum.size()*(Sm_col_offset+MPO_subcol)+col,x);
	//first block col
	M.entry(modelvertex.Spectrum.size()*(Sm_row_offset+MPO_subrow)+i,col,x*(halfJxy+momentumfactor));
	//lowest block row
	M.entry(offset_to_last_block+col,modelvertex.Spectrum.size()*(Sp_col_offset+MPO_subcolSp)+i,conj(x));
	//first block col
	M.entry(modelvertex.Spectrum.size()*(Sp_row_offset+MPO_subrowSp)+col,i,conj(x)*(halfJxy+momentumfactor));
      }
    }

    //repeat for Sz
    ajaj::Sparseint Sz_col_offset=2*differencecombinations.size()+1; //+1 for identity matrix in first block
    ajaj::Sparseint Sz_row_offset=Sz_col_offset;
    for (ajaj::Sparseint col=0;col<modelvertex.Operators[2].MatrixElements.cols();++col){
      for (ajaj::Sparseint p=modelvertex.Operators[2].MatrixElements.get_p(col);p<modelvertex.Operators[2].MatrixElements.get_p(col+1);++p){
	ajaj::Sparseint i=modelvertex.Operators[2].MatrixElements.get_i(p);
	ajaj::State diffstate=modelvertex.Spectrum[i]-modelvertex.Spectrum[col];
	//now look through charge groups to find MPO indices
	ajaj::Sparseint MPO_subcol=0;
	ajaj::Sparseint MPO_subrow=0;
	if (diffstate==-diffstate){ //involution block
	  for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	    if (diffstate==differencecombinations.InvolutionPairs[l].PairState){
	      MPO_subcol=l;
	      MPO_subrow=l; //involution pairs don't get reversed
	      break;
	    }
	  }
	}
	else {
	  for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	    if (diffstate==differencecombinations.OrderedPairs[l].PairState){
	      MPO_subcol=differencecombinations.InvolutionPairs.size()+l;
	      MPO_subrow=differencecombinations.size()-l-1;
	      break;
	    }
	  }
	}
	std::complex<double> x=modelvertex.Operators[2].MatrixElements.get_x(p);
	//lowest block row
	M.entry(offset_to_last_block+i,modelvertex.Spectrum.size()*(Sz_col_offset+MPO_subcol)+col,x);
	//first block col
	M.entry(modelvertex.Spectrum.size()*(Sz_row_offset+MPO_subrow)+i,col,x*(Jz+2.0*J2z*cos(2.0*M_PI*diffstate[1]/Nxxx)));
      }
    }
    //done forming array, now do indices
    ajaj::StateArray b;
    b.push_back(ajaj::State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    for (size_t c=0;c<3;++c){ //do for each operator
      for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	b.push_back(differencecombinations.InvolutionPairs[l].PairState);
      }
      for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	b.push_back(differencecombinations.OrderedPairs[l].PairState);
      }
    }
    b.push_back(ajaj::State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    
    std::vector<ajaj::MPXIndex> indices;
    indices.push_back(ajaj::MPXIndex(1,modelvertex.Spectrum)); //sigma primed
    indices.push_back(ajaj::MPXIndex(1,b)); //b_left
    indices.push_back(ajaj::MPXIndex(0,modelvertex.Spectrum)); //sigma
    indices.push_back(ajaj::MPXIndex(0,b)); //b_right

    return  ajaj::MPO_matrix(modelvertex.Spectrum, indices, M.finalise());

  }

  //////////////////////////////////////////////////////////////
  ajaj::Vertex VertexGenerator(const ajaj::VertexParameterArray& inputs)
  {
    //initialise vertex object
    ajaj::Vertex ModelVertex(inputs);
    if (inputs[0].Value<1){std::cout << "xxx lattice length must be at least 1!" << std::endl;exit(1);}
    const ajaj::uMPXInt xxx_chain_length(round(inputs[0].Value)); /*number of sites in xxx chain*/
    //possible issue here with integer to double in original input so overwrite
    ModelVertex.Parameters[0]=ajaj::VertexParameter("L",xxx_chain_length);
    ModelVertex.ChargeRules.push_back(0); //for xxx S_Z is a Z quantum number (i.e. the integers)
    ModelVertex.ChargeRules.push_back(xxx_chain_length); //the xxx_chain_length also tells us about the crystal momentum which is a Z_N quantum number
    const ajaj::QNVector& ChargeRules=ModelVertex.ChargeRules;
    std::cout << "Length of each xxx chain is " << ModelVertex.Parameters[0].Value << std::endl;

    //need to convert this to a file name
    std::string ending1(".dat");
    std::string ending2("_2sp.dat");
    std::ostringstream spectrumnamestream;
    spectrumnamestream << "./xxx_files/Heisdmrg_states_N" << xxx_chain_length <<ending1.c_str();
    //std::cout << spectrumnamestream.str() << std::endl;
    std::ostringstream MEnamestream;
    MEnamestream << "./xxx_files/Heisdmrg_MESzmp_N" << xxx_chain_length<<ending1.c_str();
    //std::cout << MEnamestream.str() << std::endl;

    std::cout << "Generating vertex spectrum" << std::endl;

    std::ifstream spectruminfile;
    spectruminfile.open(spectrumnamestream.str().c_str(),ios::in);
    if (spectruminfile.is_open()){
      
      std::cout << spectrumnamestream.str() << std::endl;
      //populate spectrum
      ajaj::MPXInt ID,SpnNum,S,Sz,P;
      double energy;
      ajaj::QNVector charges(2,0);
      //first line is description
      std::string dump;
      spectruminfile >> dump;
      while (spectruminfile >> ID >> SpnNum >> S >> Sz >> energy >> P) {
	charges[0]=Sz;
	charges[1]=P;
	ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,charges,energy));
      }
      spectruminfile.close();
      std::cout << "End generating chain spectrum, generated " << ModelVertex.Spectrum.size() << " states" << std::endl;
    }
    else {
      std::cout << "Couldn't open state file " << spectrumnamestream.str() << std::endl; exit(1);

    }
    std::cout << "Generating matrix elements" << std::endl;

    std::ifstream MEinfile;
    MEinfile.open(MEnamestream.str().c_str(),ios::in);
    if (MEinfile.is_open()){
      std::cout << MEnamestream.str() << std::endl;
      ModelVertex.Operators.push_back(ajaj::VertexOperator("S+",ModelVertex.Spectrum.size()));      
      ModelVertex.Operators.push_back(ajaj::VertexOperator("S-",ModelVertex.Spectrum.size()));
      ModelVertex.Operators.push_back(ajaj::VertexOperator("Sz",ModelVertex.Spectrum.size()));         
      std::string dump;
      MEinfile >> dump;
      ajaj::uMPXInt row,col;
      std::complex<double> SzME,SpME,SmME;
      double M_EPS(std::numeric_limits<double>::epsilon());
      while (MEinfile >> row >> col >> SzME >> SpME >> SmME) {
	if (abs(SpME)>M_EPS){
	  ModelVertex.Operators.at(0).MatrixElements.entry(row,col,ajaj::FixComplexPrecision(SpME));
	}
	if (abs(SmME)>M_EPS){
	  ModelVertex.Operators.at(1).MatrixElements.entry(row,col,ajaj::FixComplexPrecision(SmME));
	}
	if (abs(SzME)>M_EPS){
	  ModelVertex.Operators.at(2).MatrixElements.entry(row,col,ajaj::FixComplexPrecision(SzME));
	}
      }
    }
    else {
      std::cout << "Couldn't open state file " << MEnamestream.str() << std::endl; exit(1);
    }
    MEinfile.close();

    for (ajaj::VertexOperatorArray::iterator it=ModelVertex.Operators.begin();it!=ModelVertex.Operators.end();++it){
      it->MatrixElements.finalise();
      //it->print();
    }

    ModelVertex.Operators.push_back(ajaj::VertexOperator("Vertex_Hamiltonian"));
    ModelVertex.Operators.back().MatrixElements=ajaj::SparseMatrix(ModelVertex.basis().Energies());

    return ModelVertex;
  }

}
#endif
