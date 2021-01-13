/** @file common_defs.hpp
 *Some common type definitions. 
 *These are defined here to be used in multiple files.
 */
#ifndef COMMON_H
#define COMMON_H
#include <utility> //for std::pair
#include <vector> //for std::vector
#include <complex>
#include <cs.h> //for SuiteSparse_long
//#include "sparse_interface.hpp"

/**
 * @namespace ajaj
 * Namespace to prevent clashes between names used by AJAJ and possibly any user included files e.g. by RMK
 * For example the user might call a bunch of custom written functions to generate the vertex spectrum in vertex_generator.hpp and these might have names that clash with the internal definitions.
 */
namespace ajaj {
  /** Type for storing quantum numbers. Charges are likely to be short integers Z_N, but may have a sign e.g in the case of momenta-like quantum numbers, Z.*/
  typedef short int QuantumNumberInt;
  /** Type for indices, array labels etc. 64bit integers are preferable in case indices or array sizes grow huge.*/
  typedef SuiteSparse_long MPXInt; //needs to match matrices used
  typedef size_t uMPXInt; //unsigned ptr size type so 32bit on a 32bit machine and 64bit otherwise
  /** Type used for pairing indices. Lots of functions take a pair of MPX indices as arguments.*/
  typedef std::pair<MPXInt,MPXInt> MPXPair;
  typedef std::vector<MPXPair> MPXPairVector;
  
  typedef std::pair<uMPXInt,std::complex<double> > c_specifier;
  typedef std::vector<c_specifier> c_specifier_vector;
  typedef std::vector<c_specifier_vector> c_specifier_array;

  static const std::vector<MPXPair> contract00 {{MPXPair(0,0)}};
  static const std::vector<MPXPair> contract01(1,MPXPair(0,1));
  static const std::vector<MPXPair> contract10(1,MPXPair(1,0));
  static const std::vector<MPXPair> contract20(1,MPXPair(2,0));
  static const std::vector<MPXPair> contract21(1,MPXPair(2,1));
  static const std::vector<MPXPair> contract30(1,MPXPair(3,0));
  static const std::vector<MPXPair> contract11(1,MPXPair(1,1));
  static const std::vector<MPXPair> contract12(1,MPXPair(1,2));
  static const std::vector<MPXPair> contract13(1,MPXPair(1,3));
  static const std::vector<MPXPair> contract31(1,MPXPair(3,1));
  static const std::vector<MPXPair> contract32(1,MPXPair(3,2));
  static const std::vector<MPXPair> contract51(1,MPXPair(5,1));
  static const std::vector<MPXPair> contract74(1,MPXPair(7,4));

  static const std::vector<MPXPair> contract0011 {{MPXPair(0,0), MPXPair(1,1)}};
  static const std::vector<MPXPair> contract0013 {{MPXPair(0,0), MPXPair(1,3)}};
  static const std::vector<MPXPair> contract0015 {{MPXPair(0,0), MPXPair(1,5)}};
  static const std::vector<MPXPair> contract0022 {{MPXPair(0,0), MPXPair(2,2)}};
  static const std::vector<MPXPair> contract0110 {{MPXPair(0,1), MPXPair(1,0)}};
  static const std::vector<MPXPair> contract0120 {{MPXPair(0,1), MPXPair(2,0)}};
  static const std::vector<MPXPair> contract0130 {{MPXPair(0,1), MPXPair(3,0)}};
  static const std::vector<MPXPair> contract0150 {{MPXPair(0,1), MPXPair(5,0)}};
  static const std::vector<MPXPair> contract0210 {{MPXPair(0,2), MPXPair(1,0)}};
  static const std::vector<MPXPair> contract0241 {{MPXPair(0,2), MPXPair(4,1)}};
  static const std::vector<MPXPair> contract0310 {{MPXPair(0,3), MPXPair(1,0)}};
  static const std::vector<MPXPair> contract0311 {{MPXPair(0,3), MPXPair(1,1)}};

  static const std::vector<MPXPair> contract1071 {{MPXPair(1,0), MPXPair(7,1)}};
  static const std::vector<MPXPair> contract1122 {{MPXPair(1,1), MPXPair(2,2)}};
  static const std::vector<MPXPair> contract1130 {{MPXPair(1,1), MPXPair(3,0)}};
  static const std::vector<MPXPair> contract1162 {{MPXPair(1,1), MPXPair(6,2)}};
  static const std::vector<MPXPair> contract1220 {{MPXPair(1,2), MPXPair(2,0)}};
  static const std::vector<MPXPair> contract1221 {{MPXPair(1,2), MPXPair(2,1)}};
  static const std::vector<MPXPair> contract1325 {{MPXPair(1,3), MPXPair(2,5)}};

  static const std::vector<MPXPair> contract2031 {{MPXPair(2,0), MPXPair(3,1)}};
  static const std::vector<MPXPair> contract2032 {{MPXPair(2,0), MPXPair(3,2)}};
  static const std::vector<MPXPair> contract2041 {{MPXPair(2,0), MPXPair(4,1)}};
  static const std::vector<MPXPair> contract2130 {{MPXPair(2,1), MPXPair(3,0)}};
  static const std::vector<MPXPair> contract2150 {{MPXPair(2,1), MPXPair(5,0)}};
  static const std::vector<MPXPair> contract3041 {{MPXPair(3,0), MPXPair(4,1)}};
  static const std::vector<MPXPair> contract3150 {{MPXPair(3,1), MPXPair(5,0)}};
  static const std::vector<MPXPair> contract5472 {{MPXPair(5,4), MPXPair(7,2)}};
  static const std::vector<MPXPair> contract7254 {{MPXPair(7,2), MPXPair(5,4)}};

  static const std::vector<MPXPair> contract116276 {{MPXPair(1,1), MPXPair(6,2), MPXPair(7,6)}};
  static const std::vector<MPXPair> contract116293 {{MPXPair(1,1), MPXPair(6,2), MPXPair(9,3)}};
  static const std::vector<MPXPair> contract203142 {{MPXPair(2,0), MPXPair(3,1), MPXPair(4,2)}};
  static const std::vector<MPXPair> contract001122 {{MPXPair(0,0), MPXPair(1,1), MPXPair(2,2)}};
  static const std::vector<MPXPair> contract00112233 {{MPXPair(0,0), MPXPair(1,1), MPXPair(2,2),MPXPair(3,3)}};

  static const std::vector<MPXInt> reorder10 {{1,0}};
  static const std::vector<MPXInt> reorder021 {{0,2,1}};
  static const std::vector<MPXInt> reorder102 {{1,0,2}};
  static const std::vector<MPXInt> reorder135 {{1,3,5}};
  static const std::vector<MPXInt> reorder0213 {{0,2,1,3}};
  static const std::vector<MPXInt> reorder0132 {{0,1,3,2}};
  static const std::vector<MPXInt> reorder1032 {{1,0,3,2}};
  static const std::vector<MPXInt> reorder03214 {{0,3,2,1,4}};
  static const std::vector<MPXInt> reorder12403 {{1,2,4,0,3}};
  static const std::vector<MPXInt> reorder032415 {{0,3,2,4,1,5}};
  static const std::vector<MPXInt> reorder021354 {{0,2,1,3,5,4}};
  static const std::vector<MPXInt> reorder203145 {{2,0,3,1,4,5}};

  static char LARGESTMAGNITUDE[]={'L','M','\n'};
  static char SMALLESTREAL[]={'S','R','\n'}; //lowest real part for energies

inline std::vector<double> vector_of_reals(const std::vector<std::complex<double> >& complexvec){
    std::vector<double> ans;
    ans.reserve(complexvec.size());
    for (auto& c: complexvec){
      ans.emplace_back(c.real());
    }
    return ans;
  }
  
}

#endif
