/**
 *@file TEBD_DYN_MEASURE.cpp Driver file for dynamical measurements made on TEBD files.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>
#include <vector>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "TEBD_routines.hpp"
#include "measurement.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){

  ajaj::TEBD_DYN_Args RuntimeArgs(argc,argv); //Process the command line input
  if (RuntimeArgs.is_valid()){
    ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs)); //Create the model, based on input
    ajaj::uMPXInt CHI(RuntimeArgs.chi()); //bond dimension
    double trunc(RuntimeArgs.trunc()); //allowed truncation error
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices()); //number of vertices/nodes
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order()); //set the trotter decomposition order
    ajaj::uMPXInt t2idx(RuntimeArgs.num_indx()); // ket idx is set by command line

    //TEBD_DYN_MEASURE can only measure one pair of operators per run, because different operators lead to different
    //time evolution for unequal time correlators.
    //The required operators need their MPOs built, the below is probably overkill, as it is borrowed from other drivers
    //with more measurements.
    std::vector<ajaj::NamedMPO_matrix> generated_MPOs; //storage for MPOs
    std::vector<std::string> ops({RuntimeArgs.Op1_name(), RuntimeArgs.Op2_name()});
    std::vector<ajaj::uMPXInt> operator_storage_position;//indices for the operator's position in storage
    
    for (auto&& op : ops){
      bool found(0);
      size_t index;
      for (size_t m=0;m< generated_MPOs.size();++m){//search through generated MPOs (initial size is zero, so skipped on first run)
	if (op==generated_MPOs.at(m).Name){ //if found in the list, then skip it
	  found=1;
	  index=m;
	  break;//found
	}
      }

      if (!found){ //if not already generated
	if (myModel.vertex.operator_exists(op)){ //do matrix elements already exist in vertex def?
	  std::cout << "Operator " << op << " found in vertex definition." << std::endl << std::endl;
	  found=1;
	  index=generated_MPOs.size();
	  generated_MPOs.emplace_back(op,myModel.vertex.make_one_site_operator(op));
	}
	else {
	  std::cout <<"Assuming name includes position info" <<std::endl;
	  ajaj::ShiftedOperatorInfo trial(op);
	  if (myModel.vertex.operator_exists(trial.Name)){
	    std::cout << std::endl;
	    std::cout << "Operator " << trial.Name << " found in vertex definition." << std::endl;
	    std::cout << "Requested unitary transformation using basis quantum number " << trial.WhichCharge << " and factor " << trial.Factor <<std::endl << std::endl;
	    if (trial.WhichCharge<myModel.basis().getChargeRules().size()){
	      found=1;
	      index=generated_MPOs.size();
	      generated_MPOs.emplace_back(op,myModel.vertex.make_one_site_operator(trial));
	    }
	  }
	}
      }
      //at this point if we still can't find the operator definition, we give up
      if(!found) {
	std::cout << "Operator " << op << " couldn't be found or created from predefined matrix elements." <<std::endl;
	std::cout << "Note that the name should be a name defined in your operators file, or by a built-in model." <<std::endl;
	std::cout << "It should NOT be a .SPARSEMATRIX file name!" <<std::endl;
	return 1;
      } //not found

      operator_storage_position.emplace_back(index);
      
    }

    //load the file that contains the full time list with their indices
    std::stringstream full_tfnamestream;
    full_tfnamestream << RuntimeArgs.filename() << "_Full_idxtime.dat"; // Load the file which contains the full list of indices/times
      
    typedef std::pair<size_t,double> idx_full;//typedef for convenience
    std::vector<idx_full> full_idx_times; //storage for timeslice index and time
      
    full_idx_times.emplace_back(0,0.0);
    std::ifstream full_time_file;
    full_time_file.open(full_tfnamestream.str().c_str(),ios::in);
      
    if (full_time_file.is_open()){
      std::cout << "Opened " << full_tfnamestream.str() << "." << std::endl;
      size_t idx_full;
      double time_full;
      std::string sdump1;
      full_time_file >> std::ws;
      if (full_time_file.peek()=='#'){getline(full_time_file,sdump1);}//dump the comment line
      while (full_time_file >> idx_full >> time_full){
	getline(full_time_file,sdump1); //dump the rest of the line
	full_idx_times.emplace_back(idx_full,time_full);
      }
    }
      
    //The model and operators are now loaded
    //Next load the info from the time file.
    //The name should be INPUTFILE_TEBD_NUMVERTICES_INITIALSTATENAME_Evolution.dat
    std::stringstream tfnamestream;
    tfnamestream << RuntimeArgs.filename() << "_TEBD_" << RuntimeArgs.num_vertices() << "_" << RuntimeArgs.initial_state_name() << "_Evolution.dat";
    //From this need the indices and times
    //Note that there is no way to discern what Trotter order was used!
    typedef std::pair<size_t,double> idx_dble;//typedef for convenience
    std::vector<idx_dble> idx_times; //storage for timeslice index and time
    
    idx_times.emplace_back(0,0.0);
    
    if (tfnamestream.str().empty()){
      std::cout << "Time filename empty? Something is wrong." << std::endl;
      return 1;
    }
    else {
      std::ifstream time_file;
      time_file.open(tfnamestream.str().c_str(),ios::in);
      
      if (time_file.is_open()){
	std::cout << "Opened " << tfnamestream.str() << "." << std::endl;
	size_t idx;
	double time;
	std::string sdump;
	time_file >> std::ws;
	if (time_file.peek()=='#'){getline(time_file,sdump);}//dump the comment line
	while (time_file >> idx >> time){
	  getline(time_file,sdump); //dump the rest of the line
	  idx_times.emplace_back(idx,time);
	}

	std::cout << "There are " << idx_times.size() << " (index,time) pairs." <<std::endl;
	std::cout<<"There are "<< full_idx_times.size()<< " (index,time) pairs in the full file"<<std::endl;
	std::cout << "Note, index 0 at time 0.0 is created by default. " <<std::endl << std::endl;
	
      }
      else {
	std::cout <<"Could not open specified time data file '" << tfnamestream.str() << "'." <<std::endl;
	std::cout <<"File may be missing or you needed to specify an initial state name other than '" << RuntimeArgs.initial_state_name() << "'" << std::endl ;

	return 1;
      }
    }


    //make a list of the timeslices at which the Hamiltonian changes (includes 0)
    // Here we want to take the H_MPOs from the full list full_idx_times
    std::vector<ajaj::uMPXInt> H_MPO_indices; 
    for (auto& t : full_idx_times){
      if (ajaj::check_for_H_MPO_file(ajaj::SAVEALLNAME,t.first)){
	H_MPO_indices.push_back(t.first);
      }
      else H_MPO_indices.push_back(H_MPO_indices.back());
    }
    std::cout << "H_MPO index at each time slice " <<std::endl;
    ajaj::uMPXInt hi=0;
    for (auto &h : H_MPO_indices){
      std::cout << full_idx_times[hi++].second << " Here is it " << h <<std::endl;
    }
    
    //Before continuing need to check for and load H_MPO for time slice 0.
    /*std::stringstream FirstHMPOnamestream;
      FirstHMPOnamestream << ajaj::SAVEALLNAME << "_0_H.MPO_matrix";
      ajaj::MPO_matrix H_MPO(load_MPO_matrix(FirstHMPOnamestream.str(),myModel.basis()));
      if (!H_MPO.isConsistent()){
      std::cout << "Problem with " << FirstHMPOnamestream.str() << std::endl;
      return 1;
      }
      else {
      std::cout << "Successfully loaded " << FirstHMPOnamestream.str() << std::endl;
      H_MPO.print_indices();
      }*/
    

    //set up results file
    std::ostringstream commentlinestream;
    commentlinestream << "Index, t1, t2";
    
    for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1) {
      commentlinestream << ", Re(<" << RuntimeArgs.Op1_name() <<"(" << y1 << ",t1)" << RuntimeArgs.Op2_name() <<"(" << RuntimeArgs.y2() << ",t2)>)";
      commentlinestream << ", Im(<" << RuntimeArgs.Op1_name() <<"(" << y1 << ",t1)" << RuntimeArgs.Op2_name() <<"(" << RuntimeArgs.y2() << ",t2)>)";
    }
    
    std::stringstream outfilenamestream;
    outfilenamestream<<RuntimeArgs.filename()<<"_DYN_"<<number_of_vertices<<"_"<<RuntimeArgs.initial_state_name()<<"_" << RuntimeArgs.Op1_name() << "_y1_" << RuntimeArgs.Op2_name() << "_"<< RuntimeArgs.y2() <<".dat";
    ajaj::DataOutput results(outfilenamestream.str(),commentlinestream.str());
    ajaj::DataOutput dummyresults;
    ajaj::DataOutput dummyfullresults;

    std::string WorkingName("TempMPS");
    

    ajaj::uMPXInt H_MPO_idx=0;
    ajaj::uMPXInt results_index=0;

    //check t2 index is in the list of saved states, and if so, what the index is
    ajaj::uMPXInt ketsliceidx=0;
    while(ketsliceidx<idx_times.size()){
      if (idx_times[ketsliceidx].first==t2idx) break;
      else ketsliceidx++;
    }
    
    
    if (ketsliceidx<idx_times.size()){ //only do this if requested t2 was in the list
      double t2=idx_times[ketsliceidx].second;
      std::cout << "ket idx " << t2idx << ", here t2 = " << t2 <<std::endl;

      //We will need to load the MPS at each timeslice
      //Its files should be named SAVEALLNAME_timesliceidx_INITIALSTATENAME_Left_i.MPS_matrix
      //SAVEALLNAME is defined in TEBD_routines.hpp
      
      std::stringstream mpsrootnamestream;
      mpsrootnamestream << ajaj::SAVEALLNAME << "_" << t2idx << "_" << RuntimeArgs.initial_state_name();
      
      ajaj::FiniteMPS F(myModel.basis(),mpsrootnamestream.str(),WorkingName,number_of_vertices,1/*should be canonical*/,number_of_vertices); //load finite MPS from storage, and save a working copy 'F'.

      //equal time bit for testing
      // {
      // 	ajaj::ConstFiniteMPS CF(F);
      // 	std::vector<std::complex<double> > equal_time_results;
      // 	for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1){
      // 	  equal_time_results.emplace_back(ajaj::GeneralisedOverlap(CF,std::vector<ajaj::meas_pair>{{y1,&(generated_MPOs[operator_storage_position[0]].Matrix)},{RuntimeArgs.y2(),&(generated_MPOs[operator_storage_position[1]].Matrix)}}));
      // 	}
      // 	results.push(++results_index,ajaj::Data(std::vector<double>{t2,t2},equal_time_results));
      // }
      
      //apply operator_2 to F and canonize new F
      std::complex<double> op2weight=ApplySingleVertexOperatorToMPS(generated_MPOs[operator_storage_position[1]].Matrix,F,RuntimeArgs.y2(),ajaj::MPSCanonicalType::Left);

      //need to setup a TEBD object with an initial MPO and step size

      //create TEBD object for inner loop
      std::vector<ajaj::MultiVertexMeasurement> dummy_measurements; //TEBD object requires a measurement object, even if empty.
      //time ordered part
      {
	ajaj::uMPXInt current_ket_idx=t2idx;
	ajaj::TEBD finrun(myModel.H_MPO,F,dummyresults, dummyfullresults); //doesn't mess about creating bond operators.
	
	for (ajaj::uMPXInt brasliceidx=ketsliceidx; brasliceidx<idx_times.size(); brasliceidx++){  // Loop over saved states for t1
	  ajaj::uMPXInt t1idx=idx_times[brasliceidx].first;
	  double t1=idx_times[brasliceidx].second;
	  
	  std::cout<<"bra idx is " <<t1idx<< " and t1= " << t1 << std::endl;

	  
	  //check H for this step and update if necessary
        
	  for (ajaj::uMPXInt tstepidx=current_ket_idx; tstepidx<t1idx; ++tstepidx){  //Loop over exponentials, starting from t2 and ending at t1 - 1
	    //if (t1idx == t2idx){tstepidx = t2idx;}
            std::cout<<"Index for t2 is "<<t2idx <<", index for t1 is "<<t1idx<<" and index for step is "<<tstepidx<<std::endl;

            if (tstepidx==t2idx || H_MPO_idx!=H_MPO_indices[tstepidx]){ //if first step load H, or if step when H changes load new H
                
	      std::cout<<"Here the MPO index is "<<H_MPO_indices[tstepidx] << " and the time is"<< full_idx_times[tstepidx].second <<std::endl;
            
	      std::cout<<"Delta t has a final time:" <<full_idx_times[tstepidx+1].second<<" and an initial time: "<<full_idx_times[tstepidx].second<<std::endl;
	      double stepsize = full_idx_times[tstepidx+1].second - full_idx_times[tstepidx].second ;
	      H_MPO_idx=H_MPO_indices[tstepidx];
                
	      std::cout << "Time slice " << tstepidx << ". Updating H_MPO idx to " << H_MPO_indices[tstepidx] <<std::endl;
	      std::stringstream tSliceHMPOnamestream;
	      tSliceHMPOnamestream << ajaj::SAVEALLNAME << "_" << H_MPO_idx << "_H.MPO_matrix";
	      finrun.change_bond_operator(myModel.update_H_MPO_from_file(tSliceHMPOnamestream.str()),stepsize,trotter_order);
	      //evolve one step                
            }
            finrun.evolve(1,dummy_measurements,CHI,trunc); // The real state has been evolved herec
	  }
	  current_ket_idx=t1idx;
	  
	  ajaj::ConstFiniteMPS TEBDCF(finrun.GetEvolvingState());  // We created a copy of that state to use it in order to perform unequal time measurements
        
	  std::stringstream mpst1namestream;
	  mpst1namestream << ajaj::SAVEALLNAME << "_" << t1idx << "_" << RuntimeArgs.initial_state_name(); //t1dx should refer to the coorect bra state
	  //std::cout << "Loading bra state" <<std::endl;
	  ajaj::ConstFiniteMPS BraState(myModel.basis(),mpst1namestream.str(),number_of_vertices);
	  std::vector<std::complex<double> > unequal_time_results;
	  for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1){
	    //std::cout << "measurement" <<std::endl;
	    unequal_time_results.emplace_back(op2weight*ajaj::GeneralisedOverlap(BraState,TEBDCF,generated_MPOs[operator_storage_position[0]].Matrix,y1));
	  }
          
	  results.push(++results_index,ajaj::Data(std::vector<double>{t1,t2},unequal_time_results));
	  
	}
      }
      //anti-time ordered
      if (RuntimeArgs.include_reverse()){
	ajaj::uMPXInt current_ket_idx=t2idx;
	ajaj::TEBD finrun_rev(myModel.H_MPO,F,dummyresults,dummyfullresults); //doesn't mess about creating bond operators.
	
	for (ajaj::uMPXInt brasliceidx=ketsliceidx; brasliceidx>0; brasliceidx--){  // Loop over slices for t1
	  ajaj::uMPXInt t1idx=idx_times[brasliceidx-1].first;
	  double t1=idx_times[brasliceidx-1].second;
	  std::cout<<"anti_order t1 = "<<t1<<std::endl;
	  for (ajaj::uMPXInt tstepidx=current_ket_idx; tstepidx>t1idx; --tstepidx){  //Loop over exponentials, starting from t2 and ending at t1 - 1
            //check H for this step and update if necessary
            if (tstepidx==t2idx || H_MPO_idx!=H_MPO_indices[tstepidx-1]){
	      std::cout<<"Anti-time ordeR: Here the MPO index: "<<H_MPO_indices[tstepidx-1]<<std::endl;
	      double stepsize= full_idx_times[tstepidx-1].second-full_idx_times[tstepidx].second;
	      std::cout<<"Passed here"<<std::endl;
	      H_MPO_idx=H_MPO_indices[tstepidx-1];
                
	      std::cout << "Time slice " << t1idx << ". Updating H_MPO idx to " << H_MPO_indices[H_MPO_idx] <<std::endl;
	      std::stringstream tSliceHMPOnamestream;
	      tSliceHMPOnamestream << ajaj::SAVEALLNAME << "_" << H_MPO_idx << "_H.MPO_matrix";
	      finrun_rev.change_bond_operator(myModel.update_H_MPO_from_file(tSliceHMPOnamestream.str()),stepsize,trotter_order);
                
            }
            //evolve one step
            finrun_rev.evolve(1,dummy_measurements,CHI,trunc);
	  }
	  current_ket_idx=t1idx;
	  ajaj::ConstFiniteMPS TEBDCF(finrun_rev.GetEvolvingState());
	  std::stringstream mpst1namestream;
	  mpst1namestream << ajaj::SAVEALLNAME << "_" << t1idx << "_" << RuntimeArgs.initial_state_name();
	  //std::cout << "Loading bra state" <<std::endl;
	  ajaj::ConstFiniteMPS BraState(myModel.basis(),mpst1namestream.str(),number_of_vertices);
	  std::vector<std::complex<double> > unequal_time_results;
	  for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1){
            //std::cout << "measurement" <<std::endl;
            unequal_time_results.emplace_back(op2weight*ajaj::GeneralisedOverlap(BraState,TEBDCF,generated_MPOs[operator_storage_position[0]].Matrix,y1));
	  }
	
	  results.push(++results_index,ajaj::Data(std::vector<double>{t1,t2},unequal_time_results));
	}
      }
      
    
      return 0;
    
    }
    std::cout << "Requested t2 index not in list of saved states." <<std::endl;
    return 1;
  }

  return 1;
}
