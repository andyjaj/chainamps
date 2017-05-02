/** @file MPXIndex.hpp
 * MPXIndex is a class for tensor indices.
 * It records the charges (good quantum numbers) and 'direction' of an index for conservation purposes.
 */
#ifndef MPXINDEX_H
#define MPXINDEX_H

#include "common_defs.hpp"
#include "states.hpp"
namespace ajaj {

  //forward decs for clarity
  class MPXIndex;
  MPXIndex load_MPXIndex_binary(std::ifstream& infile,const EigenStateArray& spectrum);   /**< Load an MPXIndex from binary storage format*/

  /** A directed (ingoing or outgoing) tensor index. Contains a StateArray corresponding to a matrix index, or a reference to a Basis if it is a physical index. The directional property is necessary to sort out conservation of quantum numbers. */
  class MPXIndex{
  private:
    bool m_isInward;
    bool m_isPhysical;
    StateArray m_IndexStates; //might be empty
    const StateArray* m_IndexStatesPtr; //ptr to the states
  public:
    MPXIndex(bool inward,const StateArray& indexstates) : m_isInward(inward),m_isPhysical(0),m_IndexStates(indexstates),m_IndexStatesPtr(&m_IndexStates){}; /**< Construct a matrix index */
    MPXIndex(bool inward,StateArray&& indexstates) : m_isInward(inward),m_isPhysical(0),m_IndexStates(indexstates),m_IndexStatesPtr(&m_IndexStates){}; /**< Construct a matrix index */

    MPXIndex(bool inward,const Basis& spectrum) : m_isInward(inward),m_isPhysical(1),m_IndexStatesPtr(&spectrum){}; /**< Construct a physical index. */
    MPXIndex(const MPXIndex& other) : m_isInward(other.m_isInward),m_isPhysical(other.m_isPhysical),m_IndexStates(other.m_IndexStates),m_IndexStatesPtr(m_isPhysical ? other.m_IndexStatesPtr : &m_IndexStates){}; /**< Copy constructor */
    MPXIndex(MPXIndex&& other) noexcept : m_isInward(other.m_isInward),m_isPhysical(other.m_isPhysical),m_IndexStates(std::move(other.m_IndexStates)),m_IndexStatesPtr(m_isPhysical ? other.m_IndexStatesPtr : &m_IndexStates){};/**< Move constructor */
    MPXIndex(bool inward,const MPXIndex& other) : m_isInward(inward),m_isPhysical(other.m_isPhysical),m_IndexStates(other.m_IndexStates),m_IndexStatesPtr(m_isPhysical ? other.m_IndexStatesPtr : &m_IndexStates){}; /**< Specify direction and copy states from another index. */
    //make a trivial MPXIndex with 1 state entry
    MPXIndex(const MPXIndex& other,MPXInt i) : m_isInward(other.m_isInward),m_isPhysical(0),m_IndexStates(StateArray(1,other.at(i))),m_IndexStatesPtr(&m_IndexStates){}; /** Construct a 'dummy' index using just one State from element i of other.*/

    bool Ingoing() const {return m_isInward;} /**< Returns true if index is Ingoing*/
    bool Outgoing() const {return !Ingoing();}
    bool Physical() const {return m_isPhysical;} /**< Returns true if index is Physical */
    State at(MPXInt i) const {return m_IndexStatesPtr->at(i);} /**< Returns a copy of the i-th State in the container. */
    const State& operator[](MPXInt i) const {return (*m_IndexStatesPtr)[i];}

    MPXInt size() const {return m_IndexStatesPtr->size();} /**< How many States are in the container (what is the dimension of the index?) */
    void print() const {for (StateArray::const_iterator cit=m_IndexStatesPtr->begin();cit!=m_IndexStatesPtr->end();++cit){std::cout <<*cit << std::endl;}} /**< Print all the State objects corresponding to the index */
    bool fprint_binary(std::ofstream& outfile) const;

    // MPXIndex& operator=(MPXIndex rhs){swap(*this,rhs);return *this;}
    MPXIndex& operator=(const MPXIndex& rhs){
      m_isInward=rhs.m_isInward;
      m_isPhysical=rhs.m_isPhysical;
      m_IndexStates=rhs.m_IndexStates;//let std::vector default do the work
      m_IndexStatesPtr= m_isPhysical ? rhs.m_IndexStatesPtr : &m_IndexStates;
      return *this;}
    MPXIndex& operator=(MPXIndex&& rhs){
      m_isInward=rhs.m_isInward;
      m_isPhysical=rhs.m_isPhysical;
      m_IndexStates=std::move(rhs.m_IndexStates);//let std::vector default do the work
      m_IndexStatesPtr= m_isPhysical ? rhs.m_IndexStatesPtr : &m_IndexStates;
      return *this;
    }

    MPXIndex flip() const {return MPXIndex((!this->m_isInward),*this);} /**< Switch the direction of the index */

    friend MPXIndex combine(const MPXIndex& a,const MPXIndex& b);

    friend void swap(MPXIndex& A, MPXIndex& B);
    friend bool match(const MPXIndex& A, bool conjA,const MPXIndex& B,bool conjB); /**< Check to see if two indices match up when taking a contraction. One should be Ingoing and one should be Outgoing.*/
  };

}

#endif
