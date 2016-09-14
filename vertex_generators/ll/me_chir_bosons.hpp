void me_chir(double param)
{

  double grade[2*MaxNum_a+1], temp;
  
  int n, mop, i, l, k, r;

  /* check of hermiticity of matrix elements exp(+Phi/2\beta) */
  
  for(l=0;l<num_chir;++l) { /*modes to right of op*/
    for(k=0;k<num_chir;++k) { /*modes to left of op*/
      n = 1 + cstate_nmode[l] + cstate_nmode[k];
      mop = 1 + cstate_nmode[l];
      for(r=1;r<=cstate_nmode[l];++r)
	grade[r] = -1.*cstate_mode[l][r-1];
      grade[mop] = 0.;
      for(r=1;r<=cstate_nmode[k];++r)
	grade[r+mop] = cstate_mode[k][cstate_nmode[k]-r];
      temp = element_computation(grade,n,mop,param,1.)/sqrt(norm[l]*norm[k]);
      me_plus_chi[l][k] = temp;
    }
  }

  /* check of hermiticity of matrix elements exp(\pm Phi/2\beta) */
  
  for(l=0;l<num_chir;++l) { /*modes to right of op*/
    for(k=0;k<num_chir;++k) { /*modes to left of op*/
      n = 1 + cstate_nmode[l] + cstate_nmode[k];
      mop = 1 + cstate_nmode[l];
      for(r=1;r<=cstate_nmode[l];++r)
	grade[r] = -1.*cstate_mode[l][r-1];
      grade[mop] = 0.;
      for(r=1;r<=cstate_nmode[k];++r)
	grade[r+mop] = cstate_mode[k][cstate_nmode[k]-r];
      temp = element_computation(grade,n,mop,-1.*param,1.)/sqrt(norm[l]*norm[k]);
      me_minus_chi[l][k] = temp;
      if (fabsl(me_minus_chi[l][k]-me_plus_chi[k][l]) > .000000001)
	printf("p/m me warning: %d %d %3.8f %3.8f\n",l,k,temp,me_plus_chi[k][l]);
    }
  }

  /*computation of density op matrix elements*/
  
  /*creation part of field*/
  for(l=0;l<num_chir;++l) { /*modes to right of op*/
    for(k=0;k<num_chir;++k) { /*modes to left of op*/
      me_den_creation[l][k] = 0.;
      me_phase_creation[l][k] = 0.;
      if ((cstate_nmode[l]+1)==cstate_nmode[k]) { /*left state has one more mode than right; m.e. may be non-zero*/
	for(i=1;i<UL;++i) { /*density op adds a creation a-op to ket-state; loop through possibilities*/
	  if ((cstate_lev[l]+i)==cstate_lev[k]) {/*does modified ket-state have same level as bra-state? if yes, m.e. may be non-zero*/
	    n = 2 + cstate_nmode[l] + cstate_nmode[k];
	    mop = 2 + cstate_nmode[l];
	    for(r=1;r<=cstate_nmode[l];++r)
	      grade[r] = -1.*cstate_mode[l][r-1];
	    grade[cstate_nmode[l]+1] = ((double)-1.*i);
	    grade[mop] = 0.;
	    for(r=1;r<=cstate_nmode[k];++r)
	      grade[r+mop] = cstate_mode[k][cstate_nmode[k]-r];
	    temp = element_computation(grade,n,mop,0.,1.)/sqrt(norm[l]*norm[k]);
	    me_den_creation[l][k] += temp;
	    me_phase_creation[l][k] += (temp/((double)i));
	  }
	}
      }
      //printf("test: %d %d %3.8f\n",k,l,me_den_creation[l][k]);
    }
  }

  /*destruction part of field*/
  for(l=0;l<num_chir;++l) { /*modes to right of op; |ket>*/
    for(k=0;k<num_chir;++k) { /*modes to left of op; <bra|*/
      me_den_destruction[l][k] = 0.;
      me_phase_destruction[l][k] = 0.;
      if (cstate_nmode[l]==(cstate_nmode[k]+1)) { /*right state has one more mode than left; m.e. may be non-zero*/
	for(i=1;i<UL;++i) { /*density op adds a destruction a-op to ket-state; loop through possibilities*/
	  if ((cstate_lev[l]-i)==cstate_lev[k]) {/*does modified ket-state have same level as bra-state? if yes, m.e. may be non-zero*/
	    n = 2 + cstate_nmode[l] + cstate_nmode[k];
	    mop = 2 + cstate_nmode[l];
	    for(r=1;r<=cstate_nmode[l];++r)
	      grade[r] = -1.*cstate_mode[l][r-1];
	    grade[cstate_nmode[l]+1] = ((double)i);
	    grade[mop] = 0.;
	    for(r=1;r<=cstate_nmode[k];++r)
	      grade[r+mop] = cstate_mode[k][cstate_nmode[k]-r];
	    temp = element_computation(grade,n,mop,0.,1.)/sqrt(norm[l]*norm[k]);
	    me_den_destruction[l][k] += temp;
	    me_phase_destruction[l][k] += (temp/((double)i));
	  }
	}
      }
      //printf("%d %d %3.8f\n",k,l,me_den_destruction[l][k]);
      if (fabs(me_den_destruction[l][k]-me_den_creation[k][l])>0.00000001) {
	printf("Warning: Hermiticity problem with density op me_des_%d_%d: %3.8f me_crea_%d_%d: %3.8f\n",l,k,me_den_destruction[l][k],k,l,me_den_creation[k][l]);
      }
    }
  }
}






