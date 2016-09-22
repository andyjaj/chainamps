#ifndef THEORY_ISING1D_H
#define THEORY_ISING1D_H

#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>
#include <sstream>

#include "../ajaj_common.hpp"
#include "../sparse_interface.hpp"
namespace spin_half {
ajaj::MPO_matrix MakeHamiltonian(const ajaj::Vertex& modelvertex, const ajaj::VertexParameterArray& couplingparams){
  //Lower triangular MPO
  // I   0   0
  // JK  0   0
  // HV  K'  I  //note the charges are such that K' and K differ in that Q[K[i]]+Q[K'[i]]=0

  ajaj::MPXInt lineardim=modelvertex.Spectrum.size()*3; //the actual length of the sparse matrix needed
  ajaj::MPXInt offset_to_last_block=modelvertex.Spectrum.size()*2; //offset to get to the last row of operators
  ajaj::SparseMatrix M(lineardim,lineardim,lineardim);

  //start with the really easy bits, the Identities, I, and the vertex Hamiltonian HV
  for (size_t i=0;i<modelvertex.Spectrum.size();++i){
    M.entry(i,i,1.0);
    M.entry(i+offset_to_last_block,i,modelvertex.Spectrum[i].en);
    M.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
  }
  //need the hx part too
  if (couplingparams.size()>1 && couplingparams[1].Value!=0.0){
    M.entry(0+offset_to_last_block,1,couplingparams[1].Value);
    M.entry(1+offset_to_last_block,0,couplingparams[1].Value);
  }

  M.entry(modelvertex.Spectrum.size()+1,0,1.0);
  M.entry(modelvertex.Spectrum.size()+0,1,1.0);

  M.entry(offset_to_last_block+0,modelvertex.Spectrum.size()+1,couplingparams[0].Value);
  M.entry(offset_to_last_block+1,modelvertex.Spectrum.size()+0,couplingparams[0].Value);

  ajaj::StateArray b(3,ajaj::State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    
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
  ModelVertex.ChargeRules.push_back(1); //just a dummy, a Z_1 charge...

  const ajaj::QNVector& ChargeRules=ModelVertex.ChargeRules;
  const double hz = inputs[0].Value;

  ajaj::QuantumNumberInt gsarray[1] ={0};
  ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(gsarray,gsarray+1),+hz));
  ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(gsarray,gsarray+1),-hz));

  cout << "End generating chain spectrum, generated " << ModelVertex.Spectrum.size() << " states" << endl;
  cout << "Starting Matrix Elements" << endl;
  ModelVertex.Operators.push_back(ajaj::VertexOperator("Spin x Operator",ModelVertex.Spectrum.size()));
  ModelVertex.Operators.back().MatrixElements.entry(0,1,1.0);
  ModelVertex.Operators.back().MatrixElements.entry(1,0,1.0);
  ModelVertex.Operators.back().MatrixElements.finalise();

  ModelVertex.Operators.push_back(ajaj::VertexOperator("Spin z Operator",ModelVertex.Spectrum.size()));
  ModelVertex.Operators.back().MatrixElements.entry(0,0,1.0);
  ModelVertex.Operators.back().MatrixElements.entry(1,1,-1.0);
  ModelVertex.Operators.back().MatrixElements.finalise();

  //push back Vertex Hamiltonain
  ModelVertex.Operators.push_back(ajaj::VertexOperator("Vertex_Hamiltonian"));
  ModelVertex.Operators.back().MatrixElements=ajaj::SparseMatrix(ModelVertex.basis().Energies());

  //Operators[i]
  //0 spin x
  //1 spin z
  return ModelVertex;
}

}
#endif
