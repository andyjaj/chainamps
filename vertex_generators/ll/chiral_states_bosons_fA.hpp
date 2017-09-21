/*definitions for construction of chiral portions of states*/
const int MaxPart=300; /*maximum number of partitions allowed at a given level; used in denumerating number of chiral states*/
const int MaxChiralStates=600; /*largest allowed number of chiral states*/
const int UL=10;  /*the level up to (but not including) which we keep chiral states*/
const int LL=0;

const int MaxNum_a=20; /*the maximum number of chiral a's that we allow in constructing a state, i.e.  
			 for a state, a_{-n_1} ... a_{-n_{k_L}}\bar a_{-n_1} ... \bar a_{-n_{k_R}}|N,M>,
			 both k_L and k_R must be less than MaxNum_a; we must have MaxNum_a <= MaxLev */


const int MaxLev=20; /*the largest chiral level that we can possibly consider in forming states*/

int nstates_lev[MaxLev] = {1,2,4,7,12,19,30,45,67,97,139,195,272,373,508,684,915,1212,1597,2087}; 
/*a vector containing the number of different partitions at a given level (for levels 1 through MaxLev)*/

/*variables to describe chiral components of states as well as their matrix elements*/

int num_chir, /*the number of chiral states that we obtain*/
  cstate_lev[MaxChiralStates],
  cstate_nmode[MaxChiralStates];

double me_plus_chi[MaxChiralStates][MaxChiralStates],
  me_minus_chi[MaxChiralStates][MaxChiralStates],
  me_den_creation[MaxChiralStates][MaxChiralStates],
  me_den_destruction[MaxChiralStates][MaxChiralStates],
  me_phase_creation[MaxChiralStates][MaxChiralStates],
  me_phase_destruction[MaxChiralStates][MaxChiralStates],
  norm[MaxChiralStates],
  cstate_mode[MaxChiralStates][MaxNum_a],
  cstate_norm[MaxChiralStates];

/*routine to compute chiral matrix elements*/
double element_computation(double grade[],int nm,int m,double beta,double coeff);

/*routine to determination partitions at the necessary levels*/
int determine_partition(int k, int results[MaxPart][MaxNum_a]);

void me_chir(double param);

void chiral_states()
{

  int i, j, k, l;

  int num_modes[UL][MaxChiralStates], num_part[UL], n, mop, results[UL][MaxPart][MaxNum_a];

  double grade[2*MaxNum_a+1];
    
  /*determine partition of levels in order to enumerate chiral states; one needs to determine
   the partitions at a given level N in order to determine to how to construct all possible chiral
   states, i.e. a_{-n_1} a_{-n_2} a_{-n_3} ... a_{-n_K}|N,M>, where \sum n_i = N*/
  
  for(i=0;i<MaxPart;++i) {
    for(j=0;j<MaxNum_a;++j) {
      results[0][i][j] = 0;
    }
  }
  num_part[0] = 1;
  num_modes[0][0] = 0;
  
  for(k=1;k<UL;++k) {
    for(i=0;i<MaxPart;++i) {
      for(j=0;j<MaxNum_a;++j) {
	results[k][i][j] = 0;
      }
    }
    num_part[k] = determine_partition(k,results[k]);
    for(i=0;i<num_part[k];++i) {
      for(j=0;j<MaxNum_a;++j) {
	if (results[k][i][j] == 0) {
	  num_modes[k][i] = j;
	  break;
	}
      }
    }
  }
  
  /*construct list of chiral states*/
    
  num_chir = 0;
  for(i=0;i<UL;++i) {
    for(k=0;k<num_part[i];++k) {
      cstate_lev[num_chir] = i;
      cstate_nmode[num_chir] = num_modes[i][k];
      //printf("chir_state: %d level: %d\n",num_chir,cstate_lev[num_chir]);
      for(j=0;j<num_modes[i][k];++j) {
	cstate_mode[num_chir][j] = results[i][k][j];
	//printf("%d ",results[i][k][j]);
      }
      //printf("\n");
      num_chir++;
    }
  }
  //printf("Number of chiral states (i.e. distinct sets of a's): %d\n\n",num_chir);
  
  
  /*determine norms of chiral states*/
  
  for(l=0;l<num_chir;++l) { 
    n = 1 + cstate_nmode[l] + cstate_nmode[l];
    mop = 1 + cstate_nmode[l];
    for(i=1;i<=cstate_nmode[l];++i)
      grade[i] = -1.*cstate_mode[l][cstate_nmode[l]-i];
    grade[mop] = 0.;
    for(i=1;i<=cstate_nmode[l];++i)
      grade[i+mop] = cstate_mode[l][i-1];
    norm[l] = element_computation(grade,n,mop,0.,1.);
    //printf("%d %3.8f\n",l,norm[l]);
  }    
  
}    

/*determines partition at a given conformal grade (=level); returns them in array results*/
int determine_partition(int level, int results[MaxPart][MaxNum_a])
{
  
  int i, j, k, num, results1[MaxPart][MaxNum_a], flag, total;
  
  for(i=0;i<MaxPart;++i) {
    for(j=0;j<MaxNum_a;++j) {
      results1[i][j] = 0;
    }
  }
  
  if (level==0) {
    return(0);
  }
  
  if (level==1) {
    results[0][0] = 1;
    return(1);
  }
  
    total = 1;
    results[0][0] = level;
    for(i=level-1;i>=1;i--) {
      num = determine_partition(level-i,results1);
      for(j=0;j<num;++j) {
	flag = 1;
	for(k=0;k<MaxNum_a;++k) {
	  if (results1[j][k] > i) {
	    flag = 0;
	    break;
	  }
	}
	if (flag == 1) {
	  results[total][0] = i;
	  for(k=1;k<MaxNum_a;++k)
	    results[total][k] = results1[j][k-1];
	  total++;
	}
      }
    }
    return(total);
    
}
