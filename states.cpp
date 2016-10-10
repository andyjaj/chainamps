#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <list>

#include <algorithm>
#include "ajaj_common.hpp"
#include "states.hpp"

namespace ajaj {
  State::State(const QNVector& CHARGE_RULES) : m_CHARGE_RULES(CHARGE_RULES),values(QNVector(m_CHARGE_RULES.size(),0)) {}; //fill with zeros
  State::State(const QNVector& CHARGE_RULES, QNVector v) : m_CHARGE_RULES(CHARGE_RULES),values(v) {
    if (v.size()!= CHARGE_RULES.size()){std::cout << "Charge info doesn't match rules!" << std::endl; exit(1);}
    for (size_t i=0;i<m_CHARGE_RULES.size();++i){
      if (m_CHARGE_RULES[i]>0){
	if (values[i]<0){std::cout << "Z_" << m_CHARGE_RULES[i] << " charges can't be negative! " << values[i] << " entry " << i << std::endl;exit(1);}
	values[i]=v[i] % m_CHARGE_RULES[i];
      }
    }
  };

  bool State::is_Identity() const {
    bool flag=1; //assume true
    for (size_t i=0;i<m_CHARGE_RULES.size();++i){
      if (values[i]!=0){
	flag=0;
	break;
      }
    }
    return flag;
  }

  void State::print() const {
    for (size_t i=0;i<m_CHARGE_RULES.size();++i){
      if (m_CHARGE_RULES[i]<=0){
	std::cout << "Z: " << values[i] << " ";
      }
      else {
	std::cout << "Z_" << m_CHARGE_RULES[i] << ": " << values[i] << " ";
      }
    }
    std::cout << std::endl;
  }

  void State::fprint(std::ofstream& outfile) const {
    for (size_t i=0;i<m_CHARGE_RULES.size();++i){
      outfile << values[i] << " ";
    }
    outfile << std::endl;
  }

 void State::fprint_binary(std::ofstream& outfile) const {
   outfile.write(reinterpret_cast<const char*>(&(values[0])),sizeof(QuantumNumberInt)*m_CHARGE_RULES.size());
  }

  std::istream &operator>>(std::istream &input, State &S){
   for (size_t i=0;i<S.m_CHARGE_RULES.size();++i){
      input >> S.values[i];
      if (S.values[i]<0 && S.m_CHARGE_RULES[i]>0){std::cout << "Z_" << S.m_CHARGE_RULES[i] << " charges can't be negative! " << S.values[i] << std::endl;exit(1);}
   }
    return input;
  }

  std::ostream &operator<<(std::ostream &output, const State &S){
    output << "{ ";
    for (size_t i=0;i<S.m_CHARGE_RULES.size();++i){
      output << S.values[i] << " ";
    }
    output << "}";
    return output;
  }

  void swap(State& A, State& B){
    if (&(A.m_CHARGE_RULES)==&(B.m_CHARGE_RULES)){std::swap(A.values,B.values);}
    else {std::cout << "Trying to swap states with different charge rules!" << std::endl; exit(1);}
  }

  EigenState::EigenState(const QNVector& CHARGE_RULES, QNVector v,const double e) : State(CHARGE_RULES,v),en(e) {};
  // no implicit energy State creation
  EigenState::EigenState(State s) : State(s),en(0.0) {};
  EigenState::EigenState(State s,const double e) : State(s),en(e) {};

  void swap(EigenState& A, EigenState& B){
    if (&(A.m_CHARGE_RULES)==&(B.m_CHARGE_RULES)){
      std::swap(A.en,B.en);
      std::swap(A.values,B.values);
    }
    else {std::cout << "Trying to swap states with different charge rules!" << std::endl; exit(1);}
  }
  
  void EigenState::eprint() const {
    std::cout << "Energy: " << en << " ";
    print();
  }

  void EigenStateArray::push_back(EigenState estate){
    StateArray::push_back(State(estate));
    m_energies.push_back(estate.en);
  }

  EigenState EigenStateArray::operator[](size_t i) const{
    return EigenState(StateArray::at(i),m_energies.at(i));
  }

  void EigenStateArray::print() const {
    for (size_t i=0;i<m_energies.size();++i){
      (*this)[i].eprint();
    }
  }

  void Basis::fprint_binary(std::ofstream& outfile) const{
    //charge rules size
    size_t charge_rules_size(getChargeRules().size());
    outfile.write(reinterpret_cast<const char*>(&(charge_rules_size)),sizeof(size_t));
    outfile.write(reinterpret_cast<const char*>(&(getChargeRules())),sizeof(QuantumNumberInt)*charge_rules_size);
    size_t basis_size(size());
    outfile.write(reinterpret_cast<const char*>(&(basis_size)),sizeof(size_t));
    //states
    outfile.write(reinterpret_cast<const char*>(StateArray::data()),sizeof(QuantumNumberInt)*basis_size);
    //Energies
    outfile.write(reinterpret_cast<const char*>(Energies().data()),sizeof(double)*basis_size);
  }

  bool load_Basis_binary(std::ifstream& infile,QNVector& charge_rules, Basis& basis) {
    if (infile.is_open()){
      {
	size_t charge_rules_size(0);
	infile.read(reinterpret_cast<char*>(&charge_rules_size),sizeof(size_t));
	std::cout << charge_rules_size << std::endl;
	QNVector TempRules(charge_rules_size);
	infile.read(reinterpret_cast<char*>(TempRules.data()),sizeof(QuantumNumberInt)*charge_rules_size);
	if (charge_rules.size()==0){
	  std::cout << "Loading quantum number defs from file" << std::endl;
	  charge_rules=std::move(TempRules);
	  //new case, populate
	}
	else {
	  //check equivalence
	  if (charge_rules.size()!=charge_rules_size) return 1; //error
	  for (size_t i=0;i<charge_rules.size();++i){
	    if (charge_rules[i]!=TempRules[i]) return 1;
	  }
	}
      }
      //if we got this far, then charge_rules is fine and we can use it.

      size_t basis_size(0);
      infile.read(reinterpret_cast<char*>(&basis_size),sizeof(size_t));
      StateArray states;
      states.reserve(basis_size);
      //get quantum numbers
      QuantumNumberInt* QNbuffer=new QuantumNumberInt[basis_size];
      for (size_t b=0; b< basis_size;++b){
	infile.read(reinterpret_cast<char*>(QNbuffer),sizeof(QuantumNumberInt)*charge_rules.size());
	states.emplace_back(charge_rules,QNVector(QNbuffer,QNbuffer+charge_rules.size()));
      }
      delete[] QNbuffer;
      //get energies
      std::vector<double> energies(basis_size);
      infile.read(reinterpret_cast<char*>(energies.data()),sizeof(double)*basis_size);

      Basis TempBasis(std::move(states),std::move(energies));

      if (basis.size()==0){
	std::cout << "Loading basis from file" << std::endl;
	basis=std::move(TempBasis);
      }
      else{
	if (basis.size()!= basis_size) return 1;
	for (size_t b=0; b< basis_size;++b){
	  if (TempBasis[b]!=basis[b]) return 1;
	}
      }

      return 0;
    }
    else {
      return 1;
    }
  }

  void QNCombinations::getQNCombinations(const StateArray& sa,const bool negate=0) {
    //list is not great, but we need to erase elements all over the place
    //iterate through sa

    std::list<PairStateBlock> unsorted_ans;
    MPXInt lidx=0;
    for (StateArray::const_iterator cit1=sa.begin();cit1!=sa.end();++cit1){
      MPXInt ridx=0;
      for (StateArray::const_iterator cit2=sa.begin();cit2!=sa.end();++cit2){
	bool flag=0;
	State trialpairstate=negate ? (*cit1)-(*cit2) : (*cit1)+(*cit2);

	for (std::list<PairStateBlock>::iterator sp_it=unsorted_ans.begin();sp_it!=unsorted_ans.end();++sp_it){
	  //check to see if we already have this sector
	  if (trialpairstate==sp_it->PairState){
	    //if yes, then update flag and add product state to sector
	    flag=1;
	    sp_it->IndexPairs.push_back(std::pair<MPXInt,MPXInt>(lidx,ridx));
	    break;
	  }
	}
	if (flag==0){
	  unsorted_ans.push_back(PairStateBlock(trialpairstate,lidx,ridx));
	} 
	++ridx;
      }
      ++lidx;
    }

    //need to be careful, maybe not all states have a non empty matching block?
    //but if so...
    //what if there is no entry for the zero (identity) state?
    //check for a involution states first
    std::list<PairStateBlock>::iterator it0=unsorted_ans.begin();
    while (it0!= unsorted_ans.end()){
      if (it0->PairState==-it0->PairState){ //involution state
	InvolutionPairs.push_back(*it0); //save it to deque
	it0=unsorted_ans.erase(it0); //erase it from list
      }
      else {++it0;}
    }
    //have dealt with possible zero state
    //start again, all others must occur in pairs
    if (unsorted_ans.size()!=0){ 
      std::list<PairStateBlock>::iterator it=unsorted_ans.begin();
      OrderedPairs.push_front(*it); //put new test state at front of ans
      it=unsorted_ans.erase(it); //delete from list
      while (it!=unsorted_ans.end()){ //not at end of list?
	if (-(it->PairState)==OrderedPairs.begin()->PairState){ //is this the negative state?
	  OrderedPairs.push_back(*it); //put at end of ans
	  unsorted_ans.erase(it); //delete from list
	  it=unsorted_ans.begin(); //reset iterator
	  if (it==unsorted_ans.end()){break;}
	  OrderedPairs.push_front(*it); //new test state
	  it=unsorted_ans.erase(it); //delete from list
	}
	else { //if not the negative state, increment iterator
	  ++it;
	}
      }
    }
  }

  void QNCombinations::print() const {
    std::cout << "Involution Pairs: " << std::endl;
    for (std::vector<PairStateBlock>::const_iterator it=InvolutionPairs.begin();it!=InvolutionPairs.end();++it){
      std::cout << it->PairState << std::endl;
    }
    std::cout << "Ordered Standard Pairs: " << std::endl;
    for (std::deque<PairStateBlock>::const_iterator it=OrderedPairs.begin();it!=OrderedPairs.end();++it){
      std::cout << it->PairState << std::endl;
    }
  }

}
