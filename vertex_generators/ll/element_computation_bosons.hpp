double element_computation(double grade[2*MaxNum_a+1],int n,int m,double a,double coeff)
{
    
  int i, j, rightmode, flag, neg_mode[MaxLev+1], pos_mode[MaxLev+1], diff;
    
  double grade_res[n], coeff_res, temp, temp1, term;
    
  if (n>(2*MaxNum_a+1))
    printf("Warning: number of a's involved in matrix element is too large\n");
      
  if (coeff == 0.)
    return(0.);
    
  for(i=1;i<=(m-1);++i) {
    if (grade[i] > 0.01) {
      for(j=1;j<i;++j) {
	flag = 1;
	if (fabs(grade[i]+grade[j])<.001) {
	  flag = 0;
	  break;
	}
      }
      if (flag == 1)
	return(0.);
    }
  }
  
  /*move modes to the left of the operator to the right*/
  
  if (n > m) {
    
    /*mode commutes through op*/
    
    for(i=1;i<=n;++i) {
      grade_res[i] = grade[i];
    }
	
    grade_res[m] = grade[m+1];
    grade_res[m+1] = grade[m];
	
    temp = element_computation(grade_res,n,m+1,a,coeff);
    
    /*mode is annihilated upon commutation*/
    
    for(i=1;i<m;++i) {
      grade_res[i] = grade[i];
    }
    grade_res[m] = 0.;
    for(i=m+1;i<n;++i) {
      grade_res[i] = grade[i+1];
    }
    
    temp1 = element_computation(grade_res,n-1,m,a,-1.*a*coeff);
    
    return(temp+temp1);
  }
  
  /* all modes are to the right of the operator*/
    
  if (n == m) {
    
    /* find most rightward mode with non-negative grade */
    
    rightmode = 0;
    for(i=1;i<n;++i) {
      if (grade[i] > -.5) {
	rightmode = i;
	break;
      }
    }
    
    /* if all positive modes are to the left of all negative modes
       this routine will compute the matrix element quickly without
       further recursive calls*/
	
    flag = 1;
    for(i=rightmode+1;i<n;++i) {
      if (grade[i] < -0.01)
	flag = 0;
    }
    
    if (flag == 1) {
      for(i=1;i<=MaxLev;++i) {
	neg_mode[i] = pos_mode[i] = 0;
      }
      for(i=1;i<=n;++i) {
	if (grade[i] < -.01) {
	  neg_mode[(int)(-1.*grade[i])]++;
	}
	else 
	  pos_mode[(int)grade[i]]++;
      }
      for(i=1;i<=MaxLev;++i) {
	/*printf("modes of %d: %d %d\n",i,pos_mode[i],neg_mode[i]);*/
	if (pos_mode[i] > neg_mode[i])
	  return(0.);
      }
      term = coeff;
      for(i=1;i<=MaxLev;++i) {
	diff = neg_mode[i]-pos_mode[i];
	for(j=0;j<pos_mode[i];++j)
	  term *= (neg_mode[i]-j);
	term *= (pow((double)i,(double)pos_mode[i])*pow(a,(double)diff));
      }
      return(term);
    }
	       	
    
	
    /* all modes have negative grade */
    
    
    if (rightmode == 0) {
      coeff_res = coeff*pow(a,n-1.);
      return(coeff_res);
    }
	
    /* positive grade mode is adjacent to highest weight state*/
	
    if (rightmode == 1) {
      return(0.);
    }
    
    /*or it isn't*/
    
    else if (rightmode > 1) {
      
      for(i=1;i<=n;++i) {
	grade_res[i] = grade[i];
      }
      grade_res[rightmode] = grade[rightmode-1];
      grade_res[rightmode-1] = grade[rightmode];
      
      temp = element_computation(grade_res,n,m,a,coeff);	      
      
      temp1 = 0.;
      if (fabs(grade[rightmode]+grade[rightmode-1]) < .0001) {
	for(i=1;i<(rightmode-1);++i) {
	  grade_res[i] = grade[i];
	}
	coeff_res = coeff*grade[rightmode];
	
	for(i=rightmode+1;i<=n;++i) {
	  grade_res[i-2] = grade[i];
	}
	
	temp1 = element_computation(grade_res,n-2,m-2,a,coeff_res);	      
      }
	    
      return(temp+temp1);
      
    }
  }
  
  return(-1000.); /*code never should reach this line; the value of -1000 can be used as a debug tool*/
}











