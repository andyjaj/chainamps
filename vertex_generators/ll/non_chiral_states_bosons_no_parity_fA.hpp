int non_chir_bosons(double tpi_R, double Beta, double en_cutoff,bool restrict_n)
{
  int j, k, m, n, nmax, mmax, nmin, mmin, flag, ind, mom, mom_vac[MaxN][MaxM];
  
  double en_vac[MaxN][MaxM], en;
  
  
  /*delineate states of full bosons up to a given energy*/
  
  nmax = restrict_n ? 1 : (int) (sqrt(en_cutoff/tpi_R +1.0/12.0)/Beta)+1;
  nmin = -nmax+1;
 // nmax = 1; nmin = 0;
  mmax = 1; 
  mmin = 0;
  
  if ((2*mmax + 1) > MaxM) {
    printf("Warning: MaxM is not large enough.\n");
    exit(1);
  }
  if ((2*nmax + 1) > MaxN) {
    printf("Warning: MaxN is not large enough.\n");
    exit(1);
  }
  //printf("nmax: %d mmax: %d\n",nmax,mmax);
  
  if (en_cutoff/tpi_R+1.0/12.0 > UL){
    printf("Warning: UL is not large enough.\n");
    exit(1);
  }
  
  nstate = 0;
  num_metaZ = 0;
  for(n=nmin;n<nmax;++n) {
    for(m=mmin;m<mmax;++m) {
      en_vac[-nmin+n][-mmin+m] = tpi_R*(pow(n*Beta,2.)+pow(m*.5/Beta,2.)-1./12.);
      mom_vac[-nmin+n][-mmin+m] = n*m;
      for(k=0;k<num_chir;++k) { /*run over different left modes*/
	for(j=0;j<num_chir;++j) { /*run over different right modes*/
	  en = en_vac[-nmin+n][-mmin+m] + tpi_R*(cstate_lev[k]+cstate_lev[j]);
	  mom = mom_vac[-nmin+n][-mmin+m] + cstate_lev[k]-cstate_lev[j];

	  if (en <= en_cutoff) {
	    //printf("state: %ld %d %d %d %d %d %d %3.8f %d %3.8f %3.8f\n",nstate,n,m,k,j,cstate_lev[k],cstate_lev[j],en,mom,en_vac[-nmin+n][-mmin+m],tpi_R);
	    state_en[nstate] = en+nstate*1.0e-12;

	    state_Z[0][nstate] = n; /*center of mass momentum*/
	    state_Z[1][nstate] = m; /*U(1) charge*/
	    state_Z[2][nstate] = mom; /*momentum*/

	    state_chiral[0][nstate] = k; /*left chiral state*/
	    state_chiral[1][nstate] = j; /*right chiral state*/
	    state_chiral[2][nstate] = cstate_nmode[k]+cstate_nmode[j]; /*total number of left and right chiral modes*/

	    flag = 0;
	    for(ind=0;ind<nstate;++ind) {
	      if ((state_Z[0][ind] == state_Z[0][nstate]) && (state_Z[1][ind] == state_Z[1][nstate]) && (state_Z[2][ind] == state_Z[2][nstate])) {
		state_metaZ_ind[nstate] = state_metaZ_ind[ind];
		flag = 1;
		break;
	      }
	    }
	    if (flag==0) {
	      state_metaZ_ind[nstate] = num_metaZ;
	      num_metaZ++;
	    }
	    //printf("%d %d %d %d %d\n",nstate,n,m,mom,state_metaZ_ind[nstate]);

	    state_parity[0][nstate] = 0;
	    ++nstate;
	  }
	}
      }
    }
  }

  //printf("total number of non-chiral states: %ld\n",nstate);

  int max_n, max_m, max_mom;

  /*assign metaZ quantum number*/

  /*first find maximum n, m, and mom quantum numbers*/
  max_n = 0;
  max_m = 0;
  max_mom = 0;
  for(k=0;k<nstate;++k) {
    if (abs(state_Z[0][k])>max_n)
      max_n = abs(state_Z[0][k]);
    if (abs(state_Z[1][k])>max_m)
      max_m = abs(state_Z[1][k]);
    if (abs(state_Z[2][k])>max_mom)
      max_mom = abs(state_Z[2][k]);
  }
  max_n++;
  max_m++;
  max_mom++;
  //printf("%d %d %d\n",max_n,max_m,max_mom);
  for(k=0;k<nstate;++k) {
    state_metaZ[k] = state_Z[0][k] + state_Z[1][k]*max_n + state_Z[2][k]*max_n*max_m;
  }

  /*order states in terms of energy*/
	
  shell_en(nstate);
  for(k=0;k<nstate;++k) {
    //printf("state index: %d en: %3.8f N: %ld M: %ld mom: %ld metaZ: %ld metaZind: %ld chiral_right: %ld chiral_left %ld no. of a's: %ld parity: %ld parity_sign: %ld\n",k,state_en[k],state_Z[0][k],state_Z[1][k],state_Z[2][k],state_metaZ[k],state_metaZ_ind[k],state_chiral[0][k],state_chiral[1][k],state_chiral[2][k],state_parity[0][k],state_parity[1][k]);
  }
  return nstate;
}
	

void shell_en(int n)
{
  unsigned long i,j,k,inc;
  double v;
  int v_Z[3], v_chiral[3], v_parity[2], v_meta, v_meta_ind;
  inc = 1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do {
    inc /= 3;
    for(i=inc+1;i<=n;i++) {
      v=state_en[i-1];
      for(k=0;k<2;++k) {
	v_parity[k]=state_parity[k][i-1];
      }
      v_meta=state_metaZ[i-1];
      v_meta_ind=state_metaZ_ind[i-1];
      for(k=0;k<3;++k) {
	v_Z[k]=state_Z[k][i-1];
	v_chiral[k]=state_chiral[k][i-1];
      }
      j=i;
      while (state_en[j-1-inc] > v) {
	state_en[j-1]=state_en[j-1-inc];
	for(k=0;k<2;++k) {
	  state_parity[k][j-1]=state_parity[k][j-1-inc];
	}
	state_metaZ[j-1]=state_metaZ[j-1-inc];
	state_metaZ_ind[j-1]=state_metaZ_ind[j-1-inc];
	for(k=0;k<3;++k) {
	  state_Z[k][j-1]=state_Z[k][j-1-inc];
	  state_chiral[k][j-1]=state_chiral[k][j-1-inc];
	}
	j -= inc;
	if (j <= inc) 
	  break;
      }
      state_en[j-1]=v;
      for(k=0;k<2;++k) {
	state_parity[k][j-1]=v_parity[k];
      }
      state_metaZ[j-1]=v_meta;
      state_metaZ_ind[j-1]=v_meta_ind;
      for(k=0;k<3;++k) {
	state_Z[k][j-1]=v_Z[k];
	state_chiral[k][j-1]=v_chiral[k];
      }
    }
  } while (inc > 1);
}


 
