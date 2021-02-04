//#include <vector>
//#include <array>
//#include <string>
//#include <sstream>
//#include <cstdlib>
//#include <chrono>
//#include <iomanip>

//#include "arpack_interface.hpp" //arpack
//#include "sparse_interface.hpp" //SparseHED
//#include "states.hpp" //EigenStateArray etc.
//#include "MPX.hpp" //MPX_matrix etc.
//#include "DMRG_routines.hpp" //MPX_matrix etc.
#include "iDMRG_routines.hpp" //MPX_matrix etc.
//#include "data.hpp"

namespace ajaj {

  void iDMRG::run(uMPXInt number_of_steps, double convergence_criterion, uMPXInt chi, double truncation){
    if (number_of_steps==0 && convergence_criterion<=0.0){std::cout << "Need a finite number of steps OR convergence criterion" << std::endl;}
    if (convergence_criterion>=1.0){std::cout << "Need a convergence criterion <1.0" << std::endl;}
    //currently two site only
    double convergence(1.0);
    uMPXInt step_number=0;
    while (step_number<number_of_steps || (convergence>convergence_criterion && convergence_criterion>0.0)){
      Data this_step(grow_two_vertex(chi,truncation));
      //CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
      convergence=this_step.Real_measurements[2];
      ++step_number;
      output_ref_.push(this_step);
    }
    set_2();//if num vertices is only 2, need to adjust so that result is ok as input for other methods.
  }

}
