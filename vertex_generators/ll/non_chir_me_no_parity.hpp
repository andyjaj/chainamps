//k is out state, l in state, i.e. <k|psi^\dagger|l>

double non_chir_psi_dagger(int k,int l)
{

  int kl, kr, ll, lr, N_k, N_l;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  N_k = state_Z[0][k];
  N_l = state_Z[0][l];

  if (N_k==(N_l+1))
    return(me_plus_chi[kl][ll]*me_minus_chi[kr][lr]);
  else
    return(0.);
}

//k is out state, l in state, i.e. <k|psi|l>

double non_chir_psi(int k,int l)
{

  int kl, kr, ll, lr, N_k, N_l;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  N_k = state_Z[0][k];
  N_l = state_Z[0][l];

  if (N_k==(N_l-1))
    return(me_minus_chi[kl][ll]*me_plus_chi[kr][lr]);
  else
    return(0.);

}

//k is out state, l in state, i.e. <k|I*\phi(x=0)|l>
double phase_me(double Beta, int k, int l)
{

  int kl, kr, ll, lr, diff_Z0;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  /*contribution of zero mode - only happens if chiral parts of states k,l are identical*/
  if ((kl==ll) && (kr==lr)) {
    diff_Z0 = (state_Z[0][l]-state_Z[0][k]);
    if (diff_Z0==0)
      return(0.0);
    else 
      return(pow(-1.0,(double)(diff_Z0))/diff_Z0);
  }
  /*contribution of oscillator modes of phase - only is non-zero if states have identical N,M quantum numbers but different sets of oscillator modes*/
  else if ((state_Z[0][k]==state_Z[0][l]) && (state_Z[1][k]==state_Z[1][l]) && ((kl!=ll) || (kr!=lr))) {
    if (kl==ll) /*chiral parts of states are the same*/
      return(Beta*(-1.0*me_phase_creation[kr][lr]+me_phase_destruction[kr][lr]));
    else if (kr==lr) /*anti-chiral parts of states are the same*/
      return(Beta*(-1.0*me_phase_creation[kl][ll]+me_phase_destruction[kl][ll]));
    else
      return(0.0);
  }
  else
    return(0.0);
}


//k is out state, l in state, i.e. <k|\rho(x=0)|l>
double density_me(double Beta, int k, int l)
{

  int kl, kr, ll, lr;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  /*check if states k,l are identical*/
  if (k==l) 
    return(2.*Beta*state_Z[0][k]);
  /*if not this, check if states have identical N,M quantum numbers but different oscillator modes*/
  else if ((state_Z[0][k]==state_Z[0][l]) && (state_Z[1][k]==state_Z[1][l])) {
    if (kl==ll) /*chiral parts of states are the same*/
      return(me_den_creation[kr][lr]+me_den_destruction[kr][lr]);
    if (kr==lr) /*anti-chiral parts of states are the same*/
      return(me_den_creation[kl][ll]+me_den_destruction[kl][ll]);
    else
      return(0.0);
  }
  else
    return(0.0);

}

//k is out state, l in state, i.e. <k|\int dx\rho(x)|l>
double density_int_me(double Beta, int k, int l)
{

  int kl, kr, ll, lr;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  /*check for momentum conservation*/
  if (state_Z[2][k]==state_Z[2][l]) { 
    /*check if states k,l are identical*/
    if (k==l)
      return(2.*Beta*state_Z[0][k]);
    /*if not this, check if states have identical N,M quantum numbers but different oscillator modes*/
    else if ((state_Z[0][k]==state_Z[0][l]) && (state_Z[1][k]==state_Z[1][l]) && ((kl!=ll) || (kr!=lr))) {
      if (kl==ll) /*chiral parts of states are the same -- but with momentum conservation this should be zero automatically*/
	return(me_den_creation[kr][lr]+me_den_destruction[kr][lr]);
      if (kr==lr) /*anti-chiral parts of states are the same -- but with momentum conservation this should be zero automatically*/
	return(me_den_creation[kl][ll]+me_den_destruction[kl][ll]);
      else
	return(0.0);
    }
    else
      return(0.0);
  }
  else {
    return(0.0);
  }

}


/*matrix element <k|a_{n}|l> for the left a_{n} mode*/

double a_n(int k,int l,int n)
{

  int i, kl, kr, ll, lr, M_k, M_l, N_k, N_l, mom_k, mom_l;

  int count_ll[UL], count_kl[UL], temp_int;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  M_k = state_Z[1][k];
  M_l = state_Z[1][l];

  N_k = state_Z[0][k];
  N_l = state_Z[0][l];

  mom_k = state_Z[2][k];
  mom_l = state_Z[2][l];

  //printf("k: %d l: %d M_k: %d M_l: %d N_k: %d N_l: %d mk: %d ml: %d kr: %d kl: %d\n",k,l,M_k,M_l,N_k,N_l,mom_k,mom_l,kr,lr);
  /*check if momentum, N, M, anti-chiral conservation is obeyed*/
  if ((M_k != M_l) || (N_k != N_l) || (mom_k != (mom_l-n)) || (kr != lr)) {
    return(0.0);
  }

  /*is a_n a destruction operator?*/

  if (n > 0) {
    if (cstate_nmode[kl] != (cstate_nmode[ll]-1)) {
      return(0.);
    }
    else {
      for(i=0;i<UL;++i) {
	count_kl[i] = 0;
	count_ll[i] = 0;
      }
      for(i=0;i<cstate_nmode[kl];++i) {
	count_kl[(int)(cstate_mode[kl][i])]++;
      }
      for(i=0;i<cstate_nmode[ll];++i) {
	count_ll[(int)(cstate_mode[ll][i])]++;
      }
      for(i=0;i<UL;++i) {
	if (i != n) {
	  if (count_kl[i]!=count_ll[i])
	    return(0.);
	}
      }
      if ((count_kl[n]+1)==count_ll[n]) 
	return(double(n*count_ll[n]*sqrt(norm[kl]/norm[ll])));
      else
	return(0.);
    }
  }

  /*is a_n a creation operator?*/

  if (n < 0) {
    if (cstate_nmode[kl] != (cstate_nmode[ll]+1)) {
      return(0.0);
    }
    else {
      for(i=0;i<UL;++i) {
	count_kl[i] = 0;
	count_ll[i] = 0;
      }
      for(i=0;i<cstate_nmode[kl];++i) {
	temp_int = cstate_mode[kl][i];
	//printf("cstate_mode[%d][%d]=%3.3f %d\n",kl,i,cstate_mode[kl][i],temp_int);
	count_kl[temp_int]++;
      }
      for(i=0;i<cstate_nmode[ll];++i) {
	//printf("cstate_nmode[%d]=%3.3f\n",ll,cstate_mode[ll][i]);
	count_ll[(int)(cstate_mode[ll][i]+.001)]++;
      }
      for(i=0;i<UL;++i) {
	//printf("count_kl: %d %d\n",i,count_kl[i]);
	//printf("count_ll: %d %d\n",i,count_ll[i]);
	if (i != (-n)) {
	  if (count_kl[i]!=count_ll[i])
	    return(0.);
	}
      }
      if ((count_kl[-n]-1)==count_ll[-n]) {
	//return(double(-1.*n*count_kl[-n]*sqrt(norm[ll]/norm[kl])));
	return(double(sqrt(norm[kl]/norm[ll]))); //should also work
      }
      else
	return(0.);
    }
  }
  
  return(0.); /*should not reach here*/
}


/*matrix element <k|a_{n}|l> for the right \bar a_{n} mode*/

double bar_a_n(int k,int l,int n)
{

  int i, kl, kr, ll, lr, M_k, M_l, N_k, N_l, mom_k, mom_l;

  int count_lr[UL], count_kr[UL];

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  M_k = state_Z[1][k];
  M_l = state_Z[1][l];

  N_k = state_Z[0][k];
  N_l = state_Z[0][l];

  mom_k = state_Z[0][k];
  mom_l = state_Z[0][l];

  /*check if momentum, N, M, chiral conservation is obeyed*/
  if ((M_k != M_l) || (N_k != N_l) || (mom_k != (mom_l+n)) || (kl != ll))
    return(0.);

  /*is a_n a destruction operator?*/

  if (n > 0) {
    if (cstate_nmode[kr] != (cstate_nmode[lr]-1)) {
      return(0.);
    }
    else {
      for(i=0;i<UL;++i) {
	count_kr[i] = 0;
	count_lr[i] = 0;
      }
      for(i=0;i<cstate_nmode[kr];++i) {
	count_kr[(int)(cstate_mode[kr][i])]++;
      }
      for(i=0;i<cstate_nmode[lr];++i) {
	count_lr[(int)(cstate_mode[lr][i])]++;
      }
      for(i=0;i<UL;++i) {
	if (i != n) {
	  if (count_kr[i]!=count_lr[i])
	    return(0.);
	}
      }
      if ((count_kr[n]+1)==count_lr[n]) 
	return(double(n*count_lr[n]*sqrt(norm[kr]/norm[lr])));
      else
	return(0.);
    }
  }

  /*is bar a_n a creation operator?*/

  if (n < 0) {
    if (cstate_nmode[kr] != (cstate_nmode[lr]+1)) {
      return(0.);
    }
    else {
      for(i=0;i<UL;++i) {
	count_kr[i] = 0;
	count_lr[i] = 0;
      }
      for(i=0;i<cstate_nmode[kr];++i) {
	count_kr[(int)(cstate_mode[kr][i])]++;
      }
      for(i=0;i<cstate_nmode[lr];++i) {
	count_lr[(int)(cstate_mode[lr][i])]++;
      }
      for(i=0;i<UL;++i) {
	if (i != (-n)) {
	  if (count_kr[i]!=count_lr[i])
	    return(0.);
	}
      }
      if ((count_kr[-n]-1)==count_lr[-n]) {
	return(double(-n*count_kr[-n]*sqrt(norm[lr]/norm[kr])));
	//return(double(sqrt(norm[kr]/norm[lr]))); should also work
      }
      else
	return(0.);
    }
  }

    return(0.); /*should not reach here*/
}


/*matrix element <k|a_{-n}a_{n}|l> for the left modes, n>0*/

double n_a(int k,int l,int n)
{

  int i, kl, kr, ll, lr, M_k, M_l, N_k, N_l, mom_k, mom_l;

  int count_ll[UL], count_kl[UL], temp_int;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  M_k = state_Z[1][k];
  M_l = state_Z[1][l];

  N_k = state_Z[0][k];
  N_l = state_Z[0][l];

  mom_k = state_Z[2][k];
  mom_l = state_Z[2][l];

  if (n <= 0) {
    printf("Error: n needs to be > 0 but n=%d\n",n);
  }

  /*if n>=UL then m.e. must be zero as there will be no a_-n's in state*/
  if (n >= UL) {
    return(0.);
  }

  //printf("k: %d l: %d M_k: %d M_l: %d N_k: %d N_l: %d mk: %d ml: %d kr: %d kl: %d\n",k,l,M_k,M_l,N_k,N_l,mom_k,mom_l,kr,lr);
  /*check if momentum, N, M, left modes, right modes equal; otherwise m.e. is zero*/
  if ((M_k != M_l) || (N_k != N_l) || (mom_k != mom_l) || (kr != lr) || (cstate_nmode[kl] != cstate_nmode[ll])) {
    return(0.0);
  }


  for(i=0;i<UL;++i) {
    count_kl[i] = 0;
    count_ll[i] = 0;
  }
  for(i=0;i<cstate_nmode[kl];++i) {
    count_kl[(int)(cstate_mode[kl][i])]++;
  }
  for(i=0;i<cstate_nmode[ll];++i) {
    count_ll[(int)(cstate_mode[ll][i])]++;
  }
  for(i=0;i<UL;++i) {
    if (count_kl[i]!=count_ll[i])
      return(0.);
  }
  return(double(n*count_ll[n]));
}


/*matrix element <k|\bar a_{-n}\bar a_{n}|l> for the right modes, n>0*/

double bar_n_a(int k,int l,int n)
{

  int i, kl, kr, ll, lr, M_k, M_l, N_k, N_l, mom_k, mom_l;

  int count_lr[UL], count_kr[UL], temp_int;

  kl = state_chiral[0][k];
  kr = state_chiral[1][k];
  ll = state_chiral[0][l];
  lr = state_chiral[1][l];

  M_k = state_Z[1][k];
  M_l = state_Z[1][l];

  N_k = state_Z[0][k];
  N_l = state_Z[0][l];

  mom_k = state_Z[2][k];
  mom_l = state_Z[2][l];

  if (n <= 0) {
    printf("Error: n needs to be > 0 but n=%d\n",n);
  }

  /*if n>=UL then m.e. must be zero as there will be no a_-n's in state*/
  if (n >= UL) {
    return(0.);
  }

  //printf("k: %d l: %d M_k: %d M_l: %d N_k: %d N_l: %d mk: %d ml: %d kr: %d kl: %d\n",k,l,M_k,M_l,N_k,N_l,mom_k,mom_l,kr,lr);
  /*check if momentum, N, M, left modes, right modes equal; otherwise m.e. is zero*/
  if ((M_k != M_l) || (N_k != N_l) || (mom_k != mom_l) || (kl != ll) || (kr != lr)) {
    return(0.0);
  }


  for(i=0;i<UL;++i) {
    count_kr[i] = 0;
    count_lr[i] = 0;
  }
  for(i=0;i<cstate_nmode[kr];++i) {
    count_kr[(int)(cstate_mode[kr][i])]++;
  }
  for(i=0;i<cstate_nmode[lr];++i) {
    count_lr[(int)(cstate_mode[lr][i])]++;
  }
  for(i=0;i<UL;++i) {
    if (count_kr[i]!=count_lr[i])
      return(0.);
  }
  return(double(n*count_lr[n]));
}
