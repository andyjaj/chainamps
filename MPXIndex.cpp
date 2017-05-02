#include <iostream>
#include <fstream>

#include "MPXIndex.hpp"

namespace ajaj {

  bool MPXIndex::fprint_binary(std::ofstream& outfile) const{
    if (outfile.is_open()){
      outfile.write(reinterpret_cast<const char*>(&m_isInward),sizeof(bool));
      outfile.write(reinterpret_cast<const char*>(&m_isPhysical),sizeof(bool));
      size_t indexsize(m_IndexStatesPtr->size());
      outfile.write(reinterpret_cast<const char*>(&indexsize),sizeof(size_t));
      if (!m_isPhysical && m_IndexStatesPtr->size()>0){	
	for (StateArray::const_iterator cit=m_IndexStates.begin();cit!=m_IndexStates.end();++cit){
	  cit->fprint_binary(outfile);
	}
      }
      return 0;
    }
    else {
      return 1;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //friends
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  void swap(MPXIndex& A, MPXIndex& B){

    std::swap(A.m_isInward,B.m_isInward);
    std::swap(A.m_isPhysical,B.m_isPhysical);
    std::swap(A.m_IndexStates,B.m_IndexStates);
    std::swap(A.m_IndexStatesPtr,B.m_IndexStatesPtr);
    /*const StateArray* oldAPtr=A.m_IndexStatesPtr;
      const StateArray* oldBPtr=B.m_IndexStatesPtr;*/

    A.m_IndexStatesPtr=A.m_isPhysical ? A.m_IndexStatesPtr : &A.m_IndexStates;
    //oops!!
    B.m_IndexStatesPtr=B.m_isPhysical ? B.m_IndexStatesPtr : &B.m_IndexStates;

  }

  bool match(const MPXIndex& A, bool conjA, const MPXIndex& B, bool conjB){
    bool flag=1;
    if (A.Physical()!=B.Physical()){flag=0;} //both must be of same type
    else if (((A.Outgoing()!=conjA )!=(B.Ingoing()!=conjB))){flag=0;} //A should be outgoing, B should be ingoing
    else if(A.m_IndexStatesPtr->size()!=B.m_IndexStatesPtr->size()){flag=0;} //sizes should match
    if (flag==0){
      std::cout <<"Indices don't match!"<<std::endl;
      std::cout << "Conjugation flags: " << conjA << " " << conjB << std::endl;
      std::cout << "Physical? " << A.Physical() << " " << B.Physical() << std::endl;
      std::cout << "Index Dir? " << A.Outgoing() << " " << B.Ingoing() << std::endl;
      std::cout << "Sizes: " << A.size() << " " << B.size() << std::endl;
    }
#ifndef NDEBUG
    else {
      for (uMPXInt i=0;i<A.m_IndexStatesPtr->size();++i){
	if (A.m_IndexStatesPtr->at(i)!=B.m_IndexStatesPtr->at(i)){
	  std::cout << "Indices don't match because charges are not equal" << std::endl;
	  std::cout << A.m_IndexStatesPtr->at(i) << std::endl;
	  std::cout << B.m_IndexStatesPtr->at(i) << std::endl;
	  exit(1);
	}
      }
    }
#endif
    return flag;
  }

  MPXIndex combine(const MPXIndex& a,const MPXIndex& b){
    StateArray S;

    const StateArray& Ref_to_a_states=*(a.m_IndexStatesPtr);
    const StateArray& Ref_to_b_states=*(b.m_IndexStatesPtr);

    if (a.Ingoing()!=b.Ingoing()){ //different directions...
      for (auto&& b_it : Ref_to_b_states){
	for (auto&& a_it : Ref_to_a_states){
	  S.emplace_back(a_it-b_it);
	}
      }
    }
    else {
      for (auto&& b_it : Ref_to_b_states){
	for (auto&& a_it : Ref_to_a_states){
	  S.emplace_back(a_it+b_it);
	}
      }
    }

    return MPXIndex(a.Ingoing(),std::move(S));
  }

  MPXIndex load_MPXIndex_binary(std::ifstream& infile,const EigenStateArray& spectrum) {
    if (infile.is_open()){
      bool Inward;
      bool Physical;
      size_t indexsize;
      infile.read(reinterpret_cast<char*>(&Inward),sizeof(bool));
      infile.read(reinterpret_cast<char*>(&Physical),sizeof(bool));
      infile.read(reinterpret_cast<char*>(&indexsize),sizeof(size_t));
      if (Physical){ return MPXIndex(Inward,spectrum);}
      else {
	//we get charge rules from the physical spectrum
	const QNVector& cr(spectrum[0].getChargeRules());
	QuantumNumberInt ncr(cr.size());
	QuantumNumberInt* buffer=new QuantumNumberInt[ncr];
	StateArray loadedstates;
	loadedstates.reserve(indexsize);
	for (size_t s=0;s<indexsize;++s){
	  infile.read(reinterpret_cast<char*>(buffer),sizeof(QuantumNumberInt)*ncr);
	  loadedstates.push_back(State(cr,QNVector(buffer,buffer+ncr)));
	}
	delete[] buffer;
	return MPXIndex(Inward,loadedstates);
      }
    }
    else {
      return MPXIndex(0,spectrum);
    }
  }

}
