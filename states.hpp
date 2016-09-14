/** @file states.hpp
 *  Defines State and EigenState.
 *  Classes to hold quantum numbers (i.e. Z_N charges) used to organise matrix operations such that good quantum numbers are conserved, and block structure can be used to reduce the size of matrix operations.
 */
#ifndef STATES_H
#define STATES_H
#include <cstdlib>
#include <vector>
#include <deque>
#include <utility>
#include <iostream>
#include "ajaj_common.hpp"

namespace ajaj {

  //Forward declarations
  class State;
  class EigenState;
  typedef std::vector<QuantumNumberInt> QNVector; /**< Container type for quantum numbers. Uses std::vector interface. */
  void swap(State& A, State& B);
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  /** State is a container for some good quantum numbers, and the rules used to combine them.
   * States can be added and subtracted, using the CHARGE_RULES to decide the correct modulo
   * arithmetic.
   * Adding State objects is equivalent to working out the resulting total charges for the tensor product of two actual objects.
   * Subtracting states is useful when one considers directed arrows in MPS contractions.
   */
  class State {
  protected:
    const QNVector& m_CHARGE_RULES; /**< A reference to the list of N's that define the Z_N charges (symmetries) the model has. A zero or negative number is taken to mean a momentum-like quantum number. */
    QNVector values; /**< Contains the actual values of the charges. */
  public:
    State(const QNVector& CHARGE_RULES); /**< Constructs a 'dummy' state with identity like properties.*/
    State(const QNVector& CHARGE_RULES, QNVector v); /**< Constructs a state with the charges given by v.*/
    State(const State& other) : m_CHARGE_RULES(other.m_CHARGE_RULES), values(other.values) {} /**< Copy constructor.*/
    State(State&& other) noexcept : m_CHARGE_RULES(std::move(other.m_CHARGE_RULES)), values(std::move(other.values)) {} /**< Move constructor.*/
    QuantumNumberInt& operator[](const QuantumNumberInt i); /**< Access to the value of charge i.*/
    State& operator=(State rhs);
    State& operator+=(const State& rhs);
    State& operator-=(const State& rhs);
    const QNVector& getChargeRules() const {return m_CHARGE_RULES;} /**< Provides read only access to CHARGE_RULES for this state*/
    State Identity() const {return State(m_CHARGE_RULES);} /**< Returns the 'identity state'. Usually a state with all charges=0, so that adding all the charges to another state doesn't change anything
.*/
    bool is_Identity() const;
    void print() const; /**< Print the values of all the charges for this state to stdout.*/
    void fprint(std::ofstream& outfile) const; /**< Print the values of all the charges for this state to outfile.*/
    void fprint_binary(std::ofstream& outfile) const;

    friend bool operator==(const State& st1, const State& st2);
    friend bool operator==(const State& st1, const EigenState& est2);
    friend State operator+(const State& st1, const EigenState& e2);
    friend State operator-(const State& st1, const EigenState& e2);
    friend State operator-(const State& st1);
    friend std::istream &operator>>(std::istream &input, State &S);
    friend std::ostream &operator<<(std::ostream &output, const State &S);
    friend void swap(State& A, State& B); /**< Swap resources between two states. Used for basic arithmetic and copy construction purposes. */
  };
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  /** Eigenstate of the local (vertex) Hamiltonian. This is a special State object that also has a defined energy eigenvalue.
   *
   */
  class EigenState : public State {
  public:
    double en; /**< The energy eigenvalue of this state with respect to the vertex Hamiltonian.*/
    EigenState(const QNVector& CHARGE_RULES, QNVector v,const double e); /**< Constructor. Constructs a state with energy e.*/
    explicit EigenState(State s); /**< No implict energy state constructor. Sets energy to 0. */
    EigenState(State s, double e); /**< Constructs from a state, with an energy */
    void eprint() const; /**< Prints energy as well as values of quantum numbers. */

    friend bool operator==(const State& st1, const EigenState& est2);
    friend void swap(EigenState& A, EigenState& B);
    friend State operator+(const State& st1, const EigenState& e2);
    friend State operator-(const State& st1, const EigenState& e2);
  };
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  typedef std::vector<State> StateArray; /**< A container for states. Useful for MPX indices. Uses std::vector interface.*/

  /** Container for EigenState. Allows for quick construction of a vector of energies. */
  class EigenStateArray : public StateArray {
  private:
    std::vector<double> m_energies; /**< The energies of the states. */
  public:
    void push_back(EigenState estate); /**< Add an EigenState to the end of the container.*/
    EigenState operator[](size_t i) const; /**< Read only access to element i of the container. */
    void print() const; /**< Print all the EigenStates */
    const std::vector<double>& Energies() const {return m_energies;}
    const QNVector& getChargeRules() const {return at(0).getChargeRules();}
  };

  /** PairStateBlock is a list of different tensor products of pairs of states that give the same total charges.*/
  struct PairStateBlock {
    State PairState; /**< The combined charges of the tensor product */
    std::vector<MPXPair> IndexPairs; /**< Indices of the states in the pair, relative to their container (e.g. some StateArray).*/
    PairStateBlock(const State& is,const std::vector<MPXPair>& ip) : PairState(is),IndexPairs(ip){}; /**< Constructor that takes an MPXPair object.*/
    PairStateBlock(const State& is,const MPXInt l,const MPXInt r) : PairState(is) { /**< Constructor that takes two integers (l,r) labelling the states' positions in their container. */
      IndexPairs.push_back(std::pair<MPXInt,MPXInt>(l,r));
    };
    MPXInt block_size() const {return IndexPairs.size();} /**< Return the number of different state combinations in this block.*/
  };

  /** Essentially a container for PairStateBlock objects. Used when generating lower triangular MPOs. */
  class QNCombinations{
  public:
    std::vector<PairStateBlock> InvolutionPairs; /**< Pairs for which the charges of the tensor product are in involution. I.e. State A = - State A. */
    std::deque<PairStateBlock> OrderedPairs; /**< Pairs for which the charges of the tensor product are not in involution. I.e. State A != - State A  */

    QNCombinations(const StateArray& sa, const bool minus){getQNCombinations(sa,minus);} /**< Constructs all the possible PairStateBlocks for a particular StateArray, by either adding the combinations, or subtracting them if minus==1 (true). */

    void getQNCombinations(const StateArray& sa, const bool minus); /**< Helper function for the constructor that finds the different combinations. */
    void print() const; /**< Print all the pair states. */
    size_t size() const {return InvolutionPairs.size()+OrderedPairs.size();} /**< Returns the number of different PairStateBlocks. */
  };


  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  inline QuantumNumberInt& State::operator[](const QuantumNumberInt i) {
    return(values[i]);
  }
  
  inline State& State::operator=(State other){
    swap(*this,other);
    return *this;
  }

  inline bool operator==(const State& st1, const State& st2){
    bool flag=1;
    if (&(st1.m_CHARGE_RULES)!=&(st2.m_CHARGE_RULES)){std::cout << "Different charge rules!" << std::endl; exit(1);}
    for (size_t i=0;i<st1.m_CHARGE_RULES.size();++i){
      if (st1.values[i] !=st2.values[i]) {flag=0; break;}
    }
    return(flag);
  }

  inline bool operator==(const State& st1, const EigenState& est2){
    bool flag=1;
    if (&(st1.m_CHARGE_RULES)!=&(est2.m_CHARGE_RULES)){std::cout << "Different charge rules!" << std::endl; exit(1);}
    for (size_t i=0;i<st1.m_CHARGE_RULES.size();++i){
      if (st1.values[i] !=est2.values[i]) {flag=0; break;}
    } 
    return(flag);
  }

  inline bool operator!=(const State& st1, const State& st2){
    return(!(st1==st2));
  }

  inline bool operator!=(const State& st1, const EigenState& est2){
    return(!(st1==est2));
  }

  ////////////////////////////////////////////////////////////////////
  //must always ensure stored form obeys charge rules
  ////////////////////////////////////////////////////////////////////
  /** Adding a State A to a State B, leads to a new State C with charges such that charge i of C = charge i of A + charge i of B, modulo the value of the i th element in the CHARGE_RULES. */
  inline State& State::operator+=(const State& rhs){
    if (&(this->m_CHARGE_RULES)!=&(rhs.m_CHARGE_RULES)){std::cout << "Different charge rules!" << std::endl; exit(1);}
    for (size_t i=0;i<m_CHARGE_RULES.size();++i){
      if (m_CHARGE_RULES[i]<=0) { //momentum like
	this->values[i]+=rhs.values[i];
      }
      else { //Z_N like
	QuantumNumberInt temp=(this->values[i]+rhs.values[i]) % m_CHARGE_RULES[i];
	this->values[i]=temp;
      }
    }
    return *this;
  }

  inline State& State::operator-=(const State& rhs){
    if (&(this->m_CHARGE_RULES)!=&(rhs.m_CHARGE_RULES)){std::cout << "Different charge rules!" << std::endl; exit(1);}
    for (size_t i=0;i<m_CHARGE_RULES.size();++i){
      if (m_CHARGE_RULES[i]<=0) { //momentum like
	this->values[i]-=rhs.values[i];
      }
      else { //Z_N like
	QuantumNumberInt temp=(this->values[i]-rhs.values[i]);
	//now correct for negative charges
	while (temp<0){
	  temp+=m_CHARGE_RULES[i];
	}
	this->values[i]=temp % m_CHARGE_RULES[i];
      }
    }
    return *this;
  }

  inline State operator+(State lhs, const State& rhs){
    lhs+=rhs;
    return lhs;
  };

  inline State operator-(State lhs, const State& rhs){
    lhs-=rhs;
    return lhs;
  };

  inline State operator-(const State& st1){
    return State(State(st1.m_CHARGE_RULES)-st1);
  };

  inline State operator+(const State& st1, const EigenState& e2){
    if (&(st1.m_CHARGE_RULES)!=&(e2.m_CHARGE_RULES)){std::cout << "Different charge rules!" << std::endl; exit(1);}
    State ans(st1.m_CHARGE_RULES);
    for (size_t i=0;i<st1.m_CHARGE_RULES.size();++i){
      if (st1.m_CHARGE_RULES[i]<=0) { //momentum like
       ans.values[i]=st1.values[i]+e2.values[i];
      }
      else { //Z_N like
	ans.values[i]=(st1.values[i]+e2.values[i]) % st1.m_CHARGE_RULES[i];
      }
    }
    return ans;
  };

  inline State operator-(const State& st1, const EigenState& e2){
    if (&(st1.m_CHARGE_RULES)!=&(e2.m_CHARGE_RULES)){std::cout << "Different charge rules!" << std::endl; exit(1);}
    State ans(st1.m_CHARGE_RULES);
    for (size_t i=0;i<st1.m_CHARGE_RULES.size();++i){
      if (st1.m_CHARGE_RULES[i]<=0) { //momentum like
       ans.values[i]=st1.values[i]-e2.values[i];
      }
      else { //Z_N like
	QuantumNumberInt temp=(st1.values[i]-e2.values[i]);
	//now correct for negative charges
	while (temp<0){
	  temp+=st1.m_CHARGE_RULES[i];
	}
	ans.values[i]=temp % st1.m_CHARGE_RULES[i];
      }
    }
    return ans;
  };

}
#endif
