#include <complex>
#include <limits>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <cassert>
#include <cs.h>

#if defined(USETBB)
#include <tbb/tbb.h>
#define TBBNZ 10000 //used as part of the threading cutoff
#endif

#include "sparse_interface.hpp"
#include "dense_interface.hpp"
#include "arpack_interface.hpp"

namespace ajaj {

  void order_rows(cs_cl* A);
  SparseMatrix&& SparseMatrix::order_rows(){ajaj::order_rows(m_array); return std::move(*this);}

  //for stripping out tiny (spurious) imaginary parts
  std::complex<double> FixComplexPrecision(std::complex<double>& C){
    if (abs(C.imag()/C.real()) < std::numeric_limits<double>::epsilon()) {return C.real();}
    else if (abs(C.real()/C.imag())< std::numeric_limits<double>::epsilon()) {return std::complex<double>(0.0,C.imag());}
    else {return C;}
  }

  inline bool singular_value_compare(const std::pair<Sparseint,double>& p1, const std::pair<Sparseint,double>& p2){return (p1.second>p2.second);}

  double Sum(const std::vector<double>& Values){
    double factor=0.0;
    for (std::vector<double>::const_iterator cit=Values.begin();cit!=Values.end();++cit){
      factor+=(*cit);
    }
    return factor;
  }

  double SquareSum(const std::vector<double>& Values){
    double factor=0.0;
    for (std::vector<double>::const_iterator cit=Values.begin();cit!=Values.end();++cit){
      factor+=(*cit) * (*cit);
    }
    return factor;
  }

  double SquareSumRescale(std::vector<double>& Values, double newsquaresum){
    double old_sq_sum=0.0;
    for (std::vector<double>::const_iterator cit=Values.begin();cit!=Values.end();++cit){
      old_sq_sum+=(*cit) * (*cit);
    }
    if (old_sq_sum<=0.0){
      std::cout << "Square sum error, value is: " << old_sq_sum<<std::endl;
      for (auto i : Values){
	std::cout << i << std::endl;
      }
      exit(1);
    }
    double factor=sqrt(newsquaresum/old_sq_sum);
    std::cout << "Square sum is: " << old_sq_sum <<", rescaling to " << newsquaresum <<std::endl;

    for (size_t i=0;i<Values.size();++i){
      Values[i]=factor*Values[i];
    }
    return old_sq_sum;
  }

  double entropy(const std::vector<double>& Values){

    // std::cout << "Sq sum = " << SquareSum(Values) << std::endl;    

    double ans=0.0;
    for (std::vector<double>::const_iterator cit=Values.begin();cit!=Values.end();++cit){
      ans-=(*cit)*(*cit)*log((*cit)*(*cit));
    }
    return ans;
  }


  void index_decompose(Sparseint idx,Sparseint numdims,const Sparseint* dimsvec,Sparseint* subidxarray){
    for (Sparseint cp=0;cp<numdims-1;++cp){
      subidxarray[cp]=idx % dimsvec[cp];
      idx/=dimsvec[cp];
    }
    subidxarray[numdims-1]=idx;
  }

  bool check_equal(const SparseMatrix& A, const SparseMatrix& B, double tol){
    if (A.rows()!=B.rows() || A.cols()!=B.cols()){ //check dimensions
      return 0;
    }
    //need to allow for spurious zeros in matrix, leading to differing numbers of 'non zeros'

    //loop through all in A, and test against B
    // if A != B then check if row/col indices match and if one of the elements is actually 0.0.

    for (Sparseint c=0;c<A.cols();++c){
      Sparseint Bp=B.get_p(c);
      Sparseint Ap=A.get_p(c);
      
      while (Ap<A.get_p(c+1) && Bp<B.get_p(c+1)){ //if we haven't reached the end of both columns
	//skip values smaller than tol
	while (abs(A.get_x(Ap))<=tol && Ap<A.get_p(c+1)){ 
	  ++Ap;
	}
	while (abs(B.get_x(Bp))<=tol && Bp<B.get_p(c+1)){ 
	  ++Bp;
	}

	if (Ap==A.get_p(c+1) && Bp==B.get_p(c+1)){
	  break; //next col
	}

	//non zeros in hand
	if (A.get_i(Ap)==B.get_i(Bp)){ //rows match?
	  if (abs(A.get_x(Ap)-B.get_x(Bp))<=tol) {
	    //equal, increment both
	    ++Ap;
	    ++Bp;
	  }
	}
	else{
	  //not equal
	  return 0; 
	}
	
	
      }
    }
      
      
    return 1;
  }


  SparseMatrix::SparseMatrix(SparseType* m_array_to_use,const bool is_it_finalised) noexcept : m_array(m_array_to_use), m_finalised(is_it_finalised) {};//constructor used by internals only

  SparseMatrix::SparseMatrix() noexcept : m_array(0),m_finalised(0) {};//useful when creating an instance that needs to be replaced later

  SparseMatrix::SparseMatrix(Sparseint param_rows,Sparseint param_cols,Sparseint param_nzmax,bool compressed) : m_array(cs_cl_spalloc(param_rows,param_cols,param_nzmax,1,0)), m_finalised(1) {
    //form directly as sparse compressed
    //i[nz]
    //x[nz]
    //p[cols+1]
  }

  SparseMatrix::SparseMatrix(Sparseint param_rows,Sparseint param_cols,Sparseint param_nzmax) : m_array(cs_cl_spalloc(param_cols,param_rows,param_nzmax,1,1)), m_finalised(0) {
    //we always form transpose initially!!!!
    //because we will take the transpose later to order the rows
    //m_array=cs_cl_spalloc(param_cols,param_rows,param_nzmax,1,1);
  }

  //default is to assume as many entries as rows
  SparseMatrix::SparseMatrix(Sparseint param_rows,Sparseint param_cols) : m_array(cs_cl_spalloc(param_cols,param_rows,param_rows,1,1)), m_finalised(0) {
    //we always form transpose initially!!!!
    //because we will take the transpose later to order the rows
    //m_array=cs_cl_spalloc(param_cols,param_rows,param_nzmax,1,1);
  }

  SparseMatrix::SparseMatrix(const SparseMatrix& other) :  m_finalised(other.m_finalised)
  {
    if (other.m_array){ //if array ptr isn't null
      //do a deep copy of the array
      //general
      bool triplet = other.m_array->nz==-1 ? 0 : 1;
      m_array=(cs_cl*)cs_cl_calloc(1,sizeof(cs_cl));
      m_array->m=(other.m_array->m);
      m_array->n=(other.m_array->n);
      m_array->nz=(other.m_array->nz);
      m_array->nzmax=(other.m_array->nzmax);
      m_array->p=(Sparseint*) cs_cl_malloc(triplet ? m_array->nzmax : m_array->n+1,sizeof(Sparseint));
      m_array->i=(Sparseint*) cs_cl_malloc(m_array->nzmax,sizeof(Sparseint));
      m_array->x = (std::complex<double>*)cs_cl_malloc(m_array->nzmax,sizeof(std::complex<double>));
      //do deep copy
      if (triplet) {      //if a triplet we need to copy nz values
	for (Sparseint i=0; i<m_array->nz;++i){
	  m_array->p[i]=other.m_array->p[i];
	  m_array->i[i]=other.m_array->i[i];
	  m_array->x[i]=other.m_array->x[i];
	}
      }
      else { //if csc then number of non zeros is other.m_array->n[other.m_array->n]
	for (Sparseint i=0; i<m_array->n+1;++i){
	  m_array->p[i]=other.m_array->p[i];
	}
	for (Sparseint i=0; i<m_array->p[m_array->n];++i){
	  m_array->i[i]=other.m_array->i[i];
	  m_array->x[i]=other.m_array->x[i];
	}
      }
      #ifndef NDEBUG
      std::cout << "Sparse COPY constructor, would be better to avoid this if possible" << std::endl;
      print_sparse_info();
      #endif
    }
    else {m_array=nullptr;
      #ifndef NDEBUG
      std::cout << "Sparse NULL COPY constructor" << std::endl;
      #endif
    }
  }

  SparseMatrix::SparseMatrix(const std::vector<complex<double> >& diag,bool inverse) :  m_array(cs_cl_spalloc(diag.size(),diag.size(),diag.size(),1,1)), m_finalised(0) //create diagonal matrix
  {
    for (size_t i=0;i<diag.size();++i){
      entry(i,i,inverse ? 1.0/diag.at(i) :diag.at(i));
    }
    cheap_finalise();
  }
  SparseMatrix::SparseMatrix(const std::vector<double>& diag,bool inverse) :  m_array(cs_cl_spalloc(diag.size(),diag.size(),diag.size(),1,1)), m_finalised(0) //create diagonal matrix
  {
    for (size_t i=0;i<diag.size();++i){
      entry(i,i,inverse ? 1.0/diag.at(i) :diag.at(i));
    }
    cheap_finalise();
  }

  SparseMatrix::~SparseMatrix(){
    cs_cl_spfree(m_array);
    m_array=nullptr;
  };


  void SparseMatrix::print_emptyness() const {if(is_finalised()) {std::cout << std::setprecision(16); std::cout << 1.0-double(get_p(cols()))/double(rows()*cols()) << std::endl;}}

  void SparseMatrix::print_l_s() const {
    if (m_finalised){
      double largest(abs(get_x(0)));
      double smallest(abs(get_x(0)));
      Sparseint p_l(0);
      Sparseint p_s(0);
      for (Sparseint p=0;p<nz();++p){
	double current(abs(get_x(p)));
	if (current > largest) {largest=current; p_l=p;}
	if (current < smallest) {smallest=current; p_s=p;}
      }
      std::cout << "Largest value: " << get_x(p_l) << " Smallest value: " << get_x(p_s) << std::endl;
    }
  }


  SparseMatrix SparseMatrix::ExtractSubMatrix(const std::pair<Sparseint,Sparseint>& RowRange, const std::pair<Sparseint,Sparseint>& ColRange) const {

    //std::cout << RowRange.first << " " << RowRange.second <<std::endl;
    //std::cout << ColRange.first << " " << ColRange.second <<std::endl;

    if (ColRange.first < 0 || RowRange.first < 0 || ColRange.second >= cols() || RowRange.second >= rows()) {
      std::cout << "Illegal range for sparse matrix extraction" <<std::endl;
      return SparseMatrix();
    }
    Sparseint num_cols(ColRange.second-ColRange.first+1);
    Sparseint num_rows(RowRange.second-RowRange.first+1);

    //std::cout << num_rows << " " << num_cols <<std::endl;
  
    //first pass, count nz
    Sparseint new_nz=0;
    for (Sparseint col=ColRange.first;col<=ColRange.second;++col){
      for (Sparseint p=get_p(col);p<get_p(col+1);++p){
	if (get_i(p) > RowRange.second) break; //beyond row range already, skip to next col
	else if (get_i(p) >= RowRange.first) ++new_nz;
	else {continue;}
      }
    }
    //std::cout << new_nz << std::endl;
    //allocate
    SparseMatrix ans(num_rows,num_cols,new_nz,1);

    //second pass, fill values, create new p, and i
    
    Sparseint count(0);
    Sparseint new_col(0);
    for (Sparseint col=ColRange.first;col<=ColRange.second;++col){
      ans.put_p(new_col)=count;
      for (Sparseint p=get_p(col);p<get_p(col+1);++p){
	if (get_i(p) > RowRange.second) break; //outside range
	else if (get_i(p) >= RowRange.first) {
	  ans.put_i(count)=get_i(p)-RowRange.first;
	  ans.put_x(count++)=get_x(p);
	}
      }
      ++new_col;
    }
    ans.put_p(new_col)=count;
    if (new_nz!=count) {std::cout << "ExtractSubMatrix error " << new_nz << " " << count << std::endl; return SparseMatrix();}
    else return ans;
  }


  SparseMatrix SparseMatrix::ExtractSubMatrix(Sparseint num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<std::pair<Sparseint,Sparseint> >& IndexVal,const bool conjugate) const {
    std::vector<std::pair<Sparseint,Sparseint> > RowIndexVal;
    std::vector<std::pair<Sparseint,Sparseint> > ColIndexVal;
    for (std::vector<std::pair<Sparseint,Sparseint> >::const_iterator cit=IndexVal.begin();cit!=IndexVal.end();++cit){
      if (cit->first < num_row_idxs){
	RowIndexVal.emplace_back(*cit);
      }
      else {
	ColIndexVal.emplace_back(*cit);
      }
    }
    //work out non trivial row and col indices
    std::vector<Sparseint> NTrowidxs(num_row_idxs);
    std::iota(NTrowidxs.begin(),NTrowidxs.end(),0);
    std::vector<Sparseint> NTcolidxs(old_idx_dims.size()-num_row_idxs);
    std::iota(NTcolidxs.begin(),NTcolidxs.end(),num_row_idxs);
    
    //work out new effective rows and cols
    for (std::vector<std::pair<Sparseint,Sparseint> >::const_iterator cit=RowIndexVal.begin();cit!=RowIndexVal.end();++cit){
      std::remove(NTrowidxs.begin(),NTrowidxs.end(),cit->first);
    }
    NTrowidxs.erase(NTrowidxs.end()-RowIndexVal.size(),NTrowidxs.end());
    for (std::vector<std::pair<Sparseint,Sparseint> >::const_iterator cit=ColIndexVal.begin();cit!=ColIndexVal.end();++cit){
      std::remove(NTcolidxs.begin(),NTcolidxs.end(),cit->first);
    }
    NTcolidxs.erase(NTcolidxs.end()-ColIndexVal.size(),NTcolidxs.end());

    Sparseint new_num_rows(1);
    Sparseint new_num_cols(1);

    for (std::vector<Sparseint>::const_iterator cit=NTrowidxs.begin();cit!=NTrowidxs.end();++cit){
      new_num_rows*=old_idx_dims[*cit];
    }
    for (std::vector<Sparseint>::const_iterator cit=NTcolidxs.begin();cit!=NTcolidxs.end();++cit){
      new_num_cols*=old_idx_dims[*cit];
    }

    SparseMatrix ans(new_num_rows,new_num_cols,nz());
    Sparseint* old_idxs=new Sparseint[old_idx_dims.size()];
    for (Sparseint old_c=0;old_c<m_array->n;++old_c){ //column
      //take apart the column index
      index_decompose(old_c,old_idx_dims.size()-num_row_idxs,old_idx_dims.data()+num_row_idxs,&(old_idxs[num_row_idxs]));
      bool colskipflag=0; //check we have a col index that matches our sub matrix
      for (std::vector<std::pair<Sparseint,Sparseint> >::const_iterator cit=ColIndexVal.begin();cit!=ColIndexVal.end();++cit){
	if (cit->second!=old_idxs[cit->first]){
	  colskipflag=1;
	  break;
	}
      }
      if (colskipflag) continue;
      else {
	for (Sparseint p=m_array->p[old_c];p<m_array->p[old_c+1];++p){ //row pointer
	  //now do the same for the row index
	  index_decompose(m_array->i[p],num_row_idxs,old_idx_dims.data(),old_idxs);
	  bool rowskipflag=0; //check we have a row index that matches our sub matrix
	  for (std::vector<std::pair<Sparseint,Sparseint> >::const_iterator cit=RowIndexVal.begin();cit!=RowIndexVal.end();++cit){
	    if (cit->second!=old_idxs[cit->first]){
	      rowskipflag=1;
	      break;
	    }
	  }
	  if (rowskipflag) continue;
	  else {
	    //ok so we have all the sub indices in hand
	    //now we have to build new row and column indices
	    Sparseint new_row=0;
	    Sparseint row_mult=1;
	    for (std::vector<Sparseint>::const_iterator cit=NTrowidxs.begin();cit!=NTrowidxs.end();++cit){
	      new_row+=row_mult*old_idxs[*cit];
	      row_mult*=old_idx_dims[*cit];
	    }   
	    Sparseint new_col=0;
	    Sparseint col_mult=1;
	    for (std::vector<Sparseint>::const_iterator cit=NTcolidxs.begin();cit!=NTcolidxs.end();++cit){
	      new_col+=col_mult*old_idxs[*cit];
	      col_mult*=old_idx_dims[*cit];
	    }
	    std::complex<double> x=m_array->x[p];
	    ans.entry(new_row,new_col,conjugate ? conj(x) : x);
	  }
	}
      }
    }
    delete[] old_idxs;
    return ans.cheap_finalise();
  }

  SparseMatrix SparseMatrix::ExtractColumns(const std::vector<Sparseint>& cols) const{
    if (m_array->nz==-1 ? 0 : 1){std::cout << "Can't extract columns from triplet, need csc" << std::endl; exit(1);}
    cs_cl* A;
    A=(cs_cl*)cs_cl_calloc(1,sizeof(cs_cl));
    A->m=(m_array->m);
    A->n=cols.size();
    A->nz=m_array->nz;
    A->p=(Sparseint*) cs_cl_malloc(cols.size()+1,sizeof(Sparseint));
    Sparseint nonzero=0;
    for (Sparseint c=0;static_cast<size_t>(c)<cols.size();++c){
      A->p[c]=nonzero;
      nonzero+=m_array->p[cols[c]+1]-m_array->p[cols[c]];
    }
    A->p[cols.size()]=nonzero;
    A->nzmax=nonzero;
    A->i=(Sparseint*) cs_cl_malloc(nonzero,sizeof(Sparseint));
    A->x = (std::complex<double>*)cs_cl_malloc(nonzero,sizeof(std::complex<double>));

    Sparseint current=0;
    for (Sparseint c=0;static_cast<size_t>(c)<cols.size();++c){
      for (Sparseint p=m_array->p[cols[c]]; p<m_array->p[cols[c]+1];++p){
	  A->i[current]=m_array->i[p];
	  A->x[current++]=m_array->x[p];
	  //++current;
      }
    }
    return SparseMatrix(A,1);
  }
  SparseMatrix SparseMatrix::ExtractColumnsAndPad(const std::vector<Sparseint>& cols,Sparseint extrarows,Sparseint extracols) const{
    if (m_array->nz==-1 ? 0 : 1){std::cout << "Can't extract columns from triplet, need csc" << std::endl; exit(1);}
    cs_cl* A;
    A=(cs_cl*)cs_cl_calloc(1,sizeof(cs_cl));
    A->m=(m_array->m)+extrarows;
    A->n=cols.size()+extracols;
    A->nz=m_array->nz;
    A->p=(Sparseint*) cs_cl_malloc(A->n+1,sizeof(Sparseint));
    Sparseint nonzero=0;
    for (Sparseint c=0;static_cast<size_t>(c)<cols.size();++c){
      A->p[c]=nonzero;
      nonzero+=m_array->p[cols[c]+1]-m_array->p[cols[c]];
    }
    for (Sparseint c=cols.size();c<=A->n;++c){
      A->p[c]=nonzero;
    }
    A->nzmax=nonzero;
    A->i=(Sparseint*) cs_cl_malloc(nonzero,sizeof(Sparseint));
    A->x = (std::complex<double>*)cs_cl_malloc(nonzero,sizeof(std::complex<double>));

    Sparseint current=0;
    for (Sparseint c=0;static_cast<size_t>(c)<cols.size();++c){
      for (Sparseint p=m_array->p[cols[c]]; p<m_array->p[cols[c]+1];++p){
	  A->i[current]=m_array->i[p];
	  A->x[current++]=m_array->x[p];
	  //++current;
      }
    }
    return SparseMatrix(A,1);
  }

  SparseMatrix SparseMatrix::ZeroLastColumns(Sparseint L) const{ //drop all finite values in the last c columns
    if (m_array->nz==-1 ? 0 : 1){std::cout << "Can't zero columns from triplet, need csc" << std::endl; exit(1);}
    if (L>=m_array->n){std::cout << "Can't zero more columns than there actually are" << std::endl; exit(1);}
    if (L<=0){return SparseMatrix(*this);}
    cs_cl* A;
    A=(cs_cl*)cs_cl_calloc(1,sizeof(cs_cl));
    A->m=(m_array->m);
    A->n=(m_array->n);
    A->nz=m_array->nz;
    A->p=(Sparseint*) cs_cl_malloc(m_array->n+1,sizeof(Sparseint));
    //copy in the first m_array->n -L columns
    for (Sparseint c=0;static_cast<size_t>(c)<m_array->n-L;++c){
      A->p[c]=m_array->p[c];
    }
    A->nzmax=m_array->p[m_array->n-L]; //nonzero = the next column ptr after the finite columns
    //now make it clear the remianing cols are empty
    for (Sparseint c=m_array->n-L;static_cast<size_t>(c)<=m_array->n;++c){
      A->p[c]=A->nzmax;
    }
    A->i=(Sparseint*) cs_cl_malloc(A->nzmax,sizeof(Sparseint));
    A->x = (std::complex<double>*)cs_cl_malloc(A->nzmax,sizeof(std::complex<double>));

    for (Sparseint c=0;static_cast<size_t>(c)<m_array->n-L;++c){
      for (Sparseint p=m_array->p[c]; p<m_array->p[c+1];++p){
	A->i[p]=m_array->i[p];
	A->x[p]=m_array->x[p];
      }
    }
    return SparseMatrix(A,1);
  }

  void SparseMatrix::loopy(void (*funcptr)(Sparseint i, Sparseint p, std::complex<double> x)){
    finalise(); //must be done first
    for (Sparseint col=0;col<m_array->n;++col){
      for (Sparseint p=m_array->p[col];p<m_array->p[col+1];++p){
        funcptr(m_array->i[p],col,m_array->x[p]);
      }
    }
  }

  void SparseMatrix::purge(){
    Sparseint i=rows();
    Sparseint j=cols();
    m_finalised=0;
    cs_cl_spfree(m_array);
    m_array=cs_cl_spalloc(j,i,10,1,1); //10 is just a small buffer size
  }

  void SparseMatrix::entry(Sparseint i, Sparseint j, std::complex<double> value){
    //if (abs(value) > SPARSETOL){
      cs_cl_entry(m_array,j,i,value); //we work with the transpose until finalised!!!!
      //}
  };

  void SparseMatrix::zero_entry(Sparseint i, Sparseint j){
    m_array->i[m_array->nz]=j;
    //row
    m_array->p[m_array->nz]=i;
    //val
    m_array->x[m_array->nz++]=0.0;
  };

  //cheap entry assumes we have preallocated correctly and have the correct dimensions!!!!
  void SparseMatrix::cheap_entry(Sparseint i, Sparseint j, std::complex<double> value){
    //col, because we entr things as transpose 
    //if (abs(value) > SPARSETOL){
      m_array->i[m_array->nz]=j;
      //row
      m_array->p[m_array->nz]=i;
      //val
      m_array->x[m_array->nz++]=value;
      //}
  };

  SparseMatrix&& SparseMatrix::finalise(){
    if (!m_finalised){
      SparseType* temp_array=m_array; 
      m_array=cs_cl_compress(temp_array); //we work with the transpose until finalised!!!!
      cs_cl_spfree(temp_array);
      cs_cl_dupl(m_array);
      //cs_cl_droptol(m_array,SPARSETOL);
      //cs_cl_dropzeros(m_array);
      temp_array=m_array;
      m_array=cs_cl_transpose(temp_array,-1);
      cs_cl_spfree(temp_array);
      m_finalised=1;
    }
    return std::move(*this);
  };

  SparseMatrix&& SparseMatrix::cheap_finalise(){ //no dupl or dropzero
    if (!m_finalised){
      SparseType* temp_array=m_array; 
      m_array=cs_cl_compress(temp_array); //we work with the transpose until finalised!!!!
      cs_cl_spfree(temp_array);
      temp_array=m_array;
      m_array=0;
      //cs_cl_droptol(temp_array,SPARSETOL);
      m_array=cs_cl_transpose(temp_array,-1);
      cs_cl_spfree(temp_array);
      m_finalised=1;
    }
    return std::move(*this);
  };

  SparseMatrix&& SparseMatrix::cheap_no_transpose_finalise(){ //compress only
    if (!m_finalised){
      SparseType* temp_array=m_array; 
      m_array=cs_cl_compress(temp_array);
      cs_cl_spfree(temp_array);
      m_finalised=1;
    }
    return std::move(*this);
  }

  SparseMatrix&& SparseMatrix::resize_finalise(Sparseint x, Sparseint y){
    if (!m_finalised){
      //remember *this is stored as its transpose
      Sparseint oldm=m_array->m;
      Sparseint oldn=m_array->n;
      m_array->m=y;
      m_array->n=x;
      SparseType* temp=cs_cl_compress(m_array);
      m_array->m=oldm;
      m_array->n=oldn;
      cs_cl_spfree(m_array);
      cs_cl_dupl(temp);
      m_array=cs_cl_transpose(temp,-1);
      cs_cl_spfree(temp);
      m_finalised=1;
    }
    return std::move(*this);
  }

  //drop elements smaller than M_EPS*largest element
  Sparseint SparseMatrix::droprelative(){
    //find largest element first
    if (!m_finalised){
      std::cout << "Array not finalised, can't drop elements!" << std::endl; exit(1);
    }
    //cs_cl_dropzeros(m_array);
    double maxabs=0.0;
    for (Sparseint elem=0;elem<m_array->p[m_array->n];++elem){
      if (abs(m_array->x[elem])>maxabs){maxabs=abs(m_array->x[elem]);}
    }
    double tol=maxabs*SPARSETOL;
    return(cs_cl_droptol(m_array,tol));
  }

  SparseMatrix& SparseMatrix::massage(){
    //find largest element first
    /*if (!m_finalised){
      std::cout << "Array not finalised, can't massage elements!" << std::endl; exit(1);
    }
    for (Sparseint a=0;a<this->nz();++a){
      if (abs(this->m_array->x[a])<10.0*SPARSETOL){this->m_array->x[a]=0.0;}
      else if (abs(imag(this->m_array->x[a]))<10.0*SPARSETOL){this->m_array->x[a]=real(this->m_array->x[a]);}
      }*/
    return *this;
  }

  SparseMatrix& SparseMatrix::transpose(){
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    SparseType* temp_array=m_array;
    m_array=cs_cl_transpose(temp_array,-1);
    cs_cl_spfree(temp_array);
    return *this;
  }

  SparseMatrix& SparseMatrix::dagger(){
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    SparseType* temp_array=m_array;
    m_array=cs_cl_transpose(temp_array,1);
    cs_cl_spfree(temp_array);
    return *this;
  }

  SparseMatrix& SparseMatrix::conjugate(){
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    for (Sparseint p=0;p<nz();++p){
      m_array->x[p]=conj(m_array->x[p]);
    }
    return *this;
  }

  SparseMatrix& SparseMatrix::rescale(std::complex<double> factor){
    //how many non zero?
    if (m_array){ //if array ptr isn't null
      //do a deep copy of the array
      //general
      //Sparseint nonzero = m_array->nz==-1 ? m_array->p[m_array->n] : m_array->nz;
      for (Sparseint i=0; i<nz();++i){
	m_array->x[i]*=factor;
      }
      return *this;
    }
    else {
      m_array=0;
      return *this;
    }
  }

  SparseMatrix SparseMatrix::copy_transpose() const{
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    //SparseType* temp_array=cs_cl_transpose(temp_array,-1);
    //return(SparseMatrix<std::complex<double> >(temp_array,1));
    return SparseMatrix(cs_cl_transpose(m_array,-1),1);
  }

  SparseMatrix SparseMatrix::copy_dagger() const{
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    //SparseType* temp_array=cs_cl_transpose(temp_array,1);
    //return(SparseMatrix<std::complex<double> >(temp_array,1));
    return SparseMatrix(cs_cl_transpose(m_array,1),1);
  }

  SparseMatrix SparseMatrix::copy_conjugate() const{
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    //SparseType* temp_array=cs_cl_transpose(temp_array,1);
    //return(SparseMatrix<std::complex<double> >(temp_array,1));
    SparseMatrix ans(copy(*this));
    for (Sparseint p=0;p<ans.nz();++p){
      ans.m_array->x[p]=conj(ans.m_array->x[p]);
    }
    return ans;
  }

  std::complex<double> SparseMatrix::trace() const{
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    if (m_array->m!=m_array->n) {std::cout << "NOT A SQUARE MATRIX" << std::endl; exit(1);}
    std::complex<double> trace=0.0;
    for (Sparseint col=0;col<m_array->n;++col){
      for (Sparseint j=m_array->p[col];j<m_array->p[col+1];++j){
	if (m_array->i[j]==col){trace+=m_array->x[j];}
	else if (m_array->i[j]>col) {break;}
      }
    }
    return trace;
  }

  std::complex<double> SparseMatrix::orthogonality_check(){
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    SparseMatrix in_dagger=copy_dagger();
    SparseMatrix check=in_dagger*(*this);
    check.droprelative();
    std::complex<double> val=check.trace()/(1.0*check.rows());
    return val;
  }

  bool SparseMatrix::is_row_ordered() const {
    if (!m_finalised) return 0;
#ifndef NDEBUG
    //check rows are ordered
    for (Sparseint c=0;c<cols();++c){
      Sparseint previous_row=-1;
      for (Sparseint p=get_p(c);p<get_p(c+1);++p){
	if (get_i(p)<=previous_row) return 0;
	else {
	  previous_row=get_i(p);
	}
      }
    }
#endif
    return 1;
  }


  std::complex<double> SparseMatrix::element(const Sparseint i, const Sparseint j) const { //super inefficient for sparse types!
    if (!m_finalised){std::cout << "Not finalised! Probably you didn't mean to do this yet!" << std::endl; exit(1);}
    for (Sparseint p=m_array->p[j];p<m_array->p[j+1];++p){
      if (m_array->i[p]==i){return m_array->x[p];}
      else if (m_array->i[p]>i){break;}
    }
    return 0.0;
  }

  std::vector<std::complex<double> > SparseMatrix::diagonal() const {
    std::vector<std::complex<double> > ans;
    for (Sparseint i=0;i<cols();++i){
      ans.push_back(element(i,i));
    }
    return ans;
  }


  SparseMatrix& SparseMatrix::operator=(SparseMatrix other){
      swap(*this,other);
      return *this;
  }

  SparseMatrix& SparseMatrix::operator*=(const SparseMatrix& rhs)//
  {
    SparseMatrix other(sparse_multiply_ajaj(*this,rhs));
    swap(*this,other);
    return *this;
    /*swap(*this,SparseMatrix(cs_cl_multiply(rhs.copy_transpose().m_array,this->copy_transpose().m_array),1).transpose());
    return *this;*/
  }
  
  SparseMatrix& SparseMatrix::operator+=(const SparseMatrix& rhs)
  {
    if(!m_finalised){std::cout << "lhs not finalised!" << std::endl;exit(1);}
    if(!rhs.m_finalised){std::cout << "rhs not finalised!" << std::endl;exit(1);}
    //have to use a temporary object because of all the pointers in m_array
    SparseMatrix temp(cs_cl_add(m_array,rhs.m_array,1.0,1.0),1);
    swap(*this,temp);
    return *this;
  }

  SparseMatrix& SparseMatrix::operator-=(const SparseMatrix& rhs)
  {
    if(!m_finalised){std::cout << "lhs not finalised!" << std::endl;exit(1);}
    if(!rhs.m_finalised){std::cout << "rhs not finalised!" << std::endl;exit(1);}
    //have to use a temporary object because of all the pointers in m_array
    SparseMatrix temp(cs_cl_add(m_array,rhs.m_array,1.0,-1.0),1);
    swap(*this,temp);
    return *this;
  }

  SparseMatrix operator*(const SparseMatrix& lhs, const SparseMatrix& rhs)//why make the copy if it is only about to be destroyed?
  {
    return sparse_multiply_ajaj(lhs,rhs);
  }

  SparseMatrix operator+(const SparseMatrix& lhs, const SparseMatrix& rhs)
  {
    return SparseMatrix(cs_cl_add(lhs.m_array,rhs.m_array,1.0,1.0),1);
  }

  SparseMatrix operator-(const SparseMatrix& lhs, const SparseMatrix& rhs)
  {
    return SparseMatrix(cs_cl_add(lhs.m_array,rhs.m_array,1.0,-1.0),1);
  }

  double SparseMatrix::norm() const {
    return cs_cl_norm(m_array);
  }

  void SparseMatrix::print() const{
    cs_cl_print(m_array,0);
  }

bool SparseMatrix::fprint(std::ofstream& outfile) const{
    if (!m_finalised){std::cout << "Not finalised! Aborting..." << std::endl;}
    if (outfile.is_open()){
      //cout << "Writing to file" << endl;
      //csi nzmax ;     /* maximum number of entries */
      //csi nz ;
      //csi m ;         /* number of rows */
      //csi n ;         /* number of columns */
      //csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
      //csi *i ;        /* row indices, size nzmax */
      //double *x ;     /* numerical values, size nzmax */

      outfile << m_array->nzmax << " " << m_array->nz << " " << m_array->m << " " << m_array->n << std::endl;

      for (Sparseint p=0;p<=m_array->n;++p){
	outfile << m_array->p[p] << std::endl;
      }
      //just output non zeros
      outfile << std::setprecision(16);
      for (Sparseint d=0;d<m_array->p[m_array->n];++d){
	outfile << m_array->i[d] << " " << m_array->x[d] << std::endl;
      }

      return(0);
    }
    else {std::cout << "SparseMatrix.fprint(), file not open" << std::endl; return(1);}
  }

  bool SparseMatrix::fprint_binary(std::ofstream& outfile) const{
    if (!m_finalised){std::cout << "Not finalised! Aborting..." << std::endl;}
    if (outfile.is_open()){
      //cout << "Writing to file" << endl;
      //csi nzmax ;     /* maximum number of entries */
      //csi nz ;
      //csi m ;         /* number of rows */
      //csi n ;         /* number of columns */
      //csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
      //csi *i ;        /* row indices, size nzmax */
      //double *x ;     /* numerical values, size nzmax */

      size_t int_size=sizeof(Sparseint);
    
      Sparseint params[4];
      params[0]=m_array->nzmax;
      params[1]=m_array->nz;
      params[2]=m_array->m;
      params[3]=m_array->n;

      outfile.write(reinterpret_cast<const char*>(params),4*int_size);
      outfile.write(reinterpret_cast<const char*>(m_array->p),(m_array->n+1)*int_size);

      if (m_array->p[m_array->n]==0){
	Sparseint i=0;
	std::complex<double> x=0.0;
	outfile.write(reinterpret_cast<const char*>(&i),int_size);
	outfile.write(reinterpret_cast<const char*>(&x),sizeof(std::complex<double>));
      }
      else{
	outfile.write(reinterpret_cast<const char*>(m_array->i),m_array->nzmax*int_size);
	outfile.write(reinterpret_cast<const char*>(m_array->x),m_array->nzmax*sizeof(std::complex<double>));
      }
      //cout << "Done writing to file" << endl;
      return(0);
    }
    else {std::cout << "SparseMatrix.fprint_binary(), binary: file not open" << std::endl; return(1);}
  }

  //friend
  SparseMatrix load_SparseMatrix_binary(std::ifstream &infile)
  {
    if (infile.is_open()){
      //cout << "Reading from file" << endl;
      //csi nzmax ;     /* maximum number of entries */
      //csi nz ;
      //csi m ;         /* number of rows */
      //csi n ;         /* number of columns */
      //csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
      //csi *i ;        /* row indices, size nzmax */
      //double *x ;     /* numerical values, size nzmax */
    
      Sparseint nzmax,nz,m,n;
      cs_cl* M;
      bool final_flag=0;

      Sparseint params[4];
      size_t int_size=sizeof(Sparseint);

      infile.read(reinterpret_cast<char*>(params),4*int_size);

      nzmax=params[0];
      nz=params[1];
      m=params[2];
      n=params[3];

      if (nz!=-1) {std::cout << "sparsefread: file not in CSC form, nz = " << nz << std::endl;exit(1);M=0;}//NULL pointer
      else{
	M=(cs_cl*)cs_cl_calloc(1,sizeof(cs_cl));
	//M=sparseallocate_empty<cs_cl>();
	M->nzmax=nzmax;
	M->nz=nz;
	M->m=m;
	M->n=n;
	if (nzmax!=0){
	  M->p=(Sparseint*)cs_dl_malloc(n+1,sizeof(Sparseint));
	  M->i=(Sparseint*)cs_dl_malloc(nzmax,sizeof(Sparseint));
	  M->x=(std::complex<double>*)cs_dl_malloc(nzmax,sizeof(std::complex<double>));

	  infile.read(reinterpret_cast<char*>(M->p),(n+1)*int_size);
	  if (M->p[n]==0){
	    M->i[0]=0;
	    M->x[0]=0.0; //dummy data to read in
	  }
	  else {
	    infile.read(reinterpret_cast<char*>(M->i),nzmax*int_size);
	    infile.read(reinterpret_cast<char*>(M->x),nzmax*sizeof(std::complex<double>));
	  }
	  final_flag=1;
	  //cout << "Done reading from file" << endl;
	  //return(M);
	}
	else { //pointers exist, even if array is empty
	  M->p=(Sparseint*)cs_dl_malloc(n+1,sizeof(Sparseint));
	  infile.read(reinterpret_cast<char*>(M->p),(n+1)*int_size);
	  M->nzmax=1; //if nzmax isn't at least 1 we tend to have problems
	  M->i=(Sparseint*)cs_dl_malloc(1,sizeof(Sparseint));
	  M->x=(std::complex<double>*)cs_dl_malloc(1,sizeof(std::complex<double>));
	  M->i[0]=0;
	  M->x[0]=0.0; //dummy data to read in
	  //cout << "Done reading from file" << endl;
	  //return(M);
	  final_flag=1;
	}
      }
      return std::move(SparseMatrix(M,final_flag));
    }
    else {std::cout << "sparsefread: file not open" << std::endl; cs_cl* M=0;return std::move(SparseMatrix(M,0));}
  }

//friend
  SparseMatrix load_SparseMatrix(std::ifstream &infile)
  {
    if (infile.is_open()){
      //cout << "Reading from file" << endl;
      //csi nzmax ;     /* maximum number of entries */
      //csi nz ;
      //csi m ;         /* number of rows */
      //csi n ;         /* number of columns */
      //csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
      //csi *i ;        /* row indices, size nzmax */
      //double *x ;     /* numerical values, size nzmax */
    
      Sparseint nzmax,nz,m,n;

      if (!(infile >> nzmax >>  nz >> m >> n)) return SparseMatrix();

      cs_cl* M=(cs_cl*)cs_cl_calloc(1,sizeof(cs_cl));
      M->nzmax=nzmax;
      M->nz=nz;
      M->m=m;
      M->n=n;
      if (nzmax<=0) nzmax=1; 

      if (nz!=-1) {
	//triplet format
	M->p=(Sparseint*)cs_dl_malloc(nzmax,sizeof(Sparseint));
	M->i=(Sparseint*)cs_dl_malloc(nzmax,sizeof(Sparseint));
	M->x=(std::complex<double>*)cs_dl_malloc(nzmax,sizeof(std::complex<double>));
	Sparseint i;
	Sparseint j;
	std::complex<double> x;
	nz=0;
	//check conditions before read in...
	while (nz<nzmax && (infile >> i >> j >> x) ){
	  //we read in as transpose, to make row ordering more efficient
	  M->p[nz]=i;
	  M->i[nz]=j;
	  M->x[nz++]=x;
	}
	if (nz!=M->nz){
	  std::cout << "ERROR: number of non zeros in file doesn't match definition! " << nz << ":" <<M->nz << std::endl;
	}
	return std::move(SparseMatrix(M,0).finalise());
      }
      else{
	M->p=(Sparseint*)cs_dl_malloc(n+1,sizeof(Sparseint));
	M->i=(Sparseint*)cs_dl_malloc(nzmax,sizeof(Sparseint));
	M->x=(std::complex<double>*)cs_dl_malloc(nzmax,sizeof(std::complex<double>));

	Sparseint col=0;
	Sparseint p;
	//check conditions before read in...
	while (col<=n && infile >>p){ //last element of p is num non zero
	  M->p[col++]=p;
	  //p first
	}

	nz=0;
	Sparseint i;
	std::complex<double> x;
	//check conditions before read in...
	while (nz < M->p[M->n] && nz<nzmax && (infile >> i >> x)){
	  M->i[nz]=i;
	  M->x[nz++]=x;
	}
	if (nz!=M->p[M->n]){
	  std::cout << "ERROR: number of non zeros in file doesn't match definition! " << nz << ":" <<M->p[M->n] << std::endl;
	  std::cout << "Returning empty" <<std::endl;
	  SparseMatrix(M,1); //cleans up M
	  return SparseMatrix();
	}

	SparseMatrix ans(M,1);

	if (!ans.is_row_ordered()){ans.order_rows();}
	
	return std::move(ans);
      }
    }
    else {std::cout << "sparsefread: file not open" << std::endl; cs_cl* M=0;return std::move(SparseMatrix(M,0));}
  }

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  SparseMatrix copy(const SparseMatrix& other){return SparseMatrix(other);}

  //is a friend function so it can check finalised status and directly access array

  SparseMatrix addmult(const SparseMatrix& A, const SparseMatrix& B,const std::complex<double> alpha,const std::complex<double> beta)
  {
    if(!A.m_finalised){std::cout << "Input array A wasn't finalised?" << std::endl; exit(1);}
    if(!B.m_finalised){std::cout << "Input array B wasn't finalised?" << std::endl; exit(1);}
    return(SparseMatrix(cs_cl_add(A.m_array,B.m_array,alpha,beta),1));
  }

  void swap (SparseMatrix& first, SparseMatrix& second) noexcept
  {
    using std::swap;
    swap(first.m_array,second.m_array);
    swap(first.m_finalised,second.m_finalised);
  }

  /*SparseMatrix reshape(const SparseMatrix& old,Sparseint newrows,Sparseint newcols)
  {
    if (!old.is_finalised()){std::cout << "Input array wasn't finalised?" << std::endl; exit(1);}
    if (old.rows()*old.cols()!=newrows*newcols){std::cout << "Incorrect reshape params" << std::endl; exit(1);}
    if (newrows==old.rows()){
      // just copy
      return SparseMatrix(old);
    }
    else {
      SparseMatrix new_matrix(newrows,newcols,old.nz());    
      for (Sparseint c=0;c<old.m_array->n;++c){
	for (Sparseint p=old.m_array->p[c];p<old.m_array->p[c+1];++p){
	  Sparseint multi=old.m_array->i[p]+(old.m_array->m)*c;
	  Sparseint new_i=multi % newrows;
	  Sparseint new_j=multi / newrows;
	  new_matrix.cheap_entry(new_i,new_j,old.m_array->x[p]);
	}
      }
      new_matrix.cheap_finalise();
      return new_matrix;
    }
  }*/

  SparseMatrix reshape(const SparseMatrix& old,Sparseint newrows)
  {
    if (!old.is_finalised()){std::cout << "Input array wasn't finalised?" << std::endl; exit(1);}
    Sparseint newcols=(old.rows()*old.cols())/newrows;
    if (old.rows()*old.cols()!=newrows*newcols){std::cout << "Incorrect reshape params" << std::endl; exit(1);}
    if (newrows==old.rows()){
      // just copy
#ifndef NDEBUG
      std::cout << "Requested reshape is the same as original, so just copying..." <<std::endl;
#endif
      return SparseMatrix(old);
    }
    else {
      //array should already be row and col ordered, so write directly into csc form
      SparseType* new_array=cs_cl_spalloc(newrows,newcols,old.nz(),1,0);
      Sparseint* numincol=new Sparseint[newcols];
      std::fill(numincol,numincol+newcols,0);
      for (Sparseint c=0;c<old.m_array->n;++c){
	for (Sparseint p=old.m_array->p[c];p<old.m_array->p[c+1];++p){
	  Sparseint multi=old.m_array->i[p]+(old.m_array->m)*c;
	  new_array->i[p]=multi % newrows;//new_i
	  numincol[multi / newrows]++; //increment col totals
	  new_array->x[p]=old.m_array->x[p]; //copy c across
	}
      }
      // now do new col pointers
      //cs_cl_cumsum(new_array->p,numincol,newcols);
      Sparseint nz=0;
      for (Sparseint z=0;z<newcols;++z){
	new_array->p[z]=nz;
	nz+=numincol[z];
      }
      new_array->p[newcols]=nz;
      delete[] numincol;
      return SparseMatrix(new_array,1);
    }
  }

  // a generalised reshape
  SparseMatrix reshape(const SparseMatrix& old,Sparseint old_num_row_idxs,Sparseint new_num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<Sparseint>& new_idx_order,const bool conjugate)
  {
    //check finalised
    if (!old.is_finalised()){std::cout << "Input array wasn't finalised?" << std::endl; exit(1);}

    //check product of dims is the same for old
    Sparseint old_prod=1;
    for (std::vector<Sparseint>::const_iterator cit=old_idx_dims.begin();cit!=old_idx_dims.end();++cit){
      old_prod*=(*cit);
    }
    if (old_prod!=old.rows()*old.cols()){std::cout << "Incorrect old params: " << old.rows()*old.cols() << ", " << old_prod << std::endl; exit(1);}
    if (old_idx_dims.size()!=new_idx_order.size()){std::cout << "Incorrect size params, new and old number of indices don't match " << old_idx_dims.size() <<" "<< new_idx_order.size() << std::endl; exit(1);}

    /*double largest(0.0);

    for (Sparseint p=0;p<old.nz();++p){
      if (abs(old.get_x(p))>largest){
	largest=abs(old.get_x(p));
      }
      }*/

    //std::cout << "LARGEST VALUE IS " << largest <<std::endl;

    Sparseint old_num_col_idxs=old_idx_dims.size()-old_num_row_idxs;
    Sparseint* multipliers=new Sparseint[old_idx_dims.size()];
    Sparseint new_num_rows=1;
    Sparseint new_num_cols=1;
    Sparseint** x1x2=0;
    Sparseint* old2new=0;
    {
      for (Sparseint i=0;i<new_num_row_idxs;++i){
	multipliers[i]=new_num_rows;
	new_num_rows*=old_idx_dims[new_idx_order[i]];
      }
      for (size_t i=new_num_row_idxs;i<new_idx_order.size();++i){
	multipliers[i]=new_num_cols;
	new_num_cols*=old_idx_dims[new_idx_order[i]];
      }
      {
	x1x2=new Sparseint*[old.rows()];
	std::fill(x1x2,x1x2+old.rows(),static_cast<Sparseint*>(0));
      }
      //make the index mulitpliers
      {
	old2new=new Sparseint[old_idx_dims.size()];
	for (size_t l=0;l<old_idx_dims.size();++l){
	  old2new[new_idx_order[l]]=l;
	}
      }
    }

    SparseMatrix new_matrix(new_num_rows,new_num_cols,old.nz());
    //create an array to store the subindices in
    Sparseint* old_col_idxs;
    Sparseint* old_row_idxs;

    old_col_idxs=new Sparseint[old_num_col_idxs];
    old_row_idxs=new Sparseint[old_num_row_idxs];
    for (Sparseint old_c=0;old_c<old.m_array->n;++old_c){ //column
      //take apart the column index if, and only if, the column is not empty
      if(old.m_array->p[old_c]!=old.m_array->p[old_c+1]){
	index_decompose(old_c,old_num_col_idxs,old_idx_dims.data()+old_num_row_idxs,old_col_idxs);
	Sparseint y1=0;
	Sparseint y2=0;
	//and now update old_idxs by their multipliers
	for (Sparseint l=0;l<old_num_col_idxs;++l){
	  if (old2new[l+old_num_row_idxs]<new_num_row_idxs){y1+=old_col_idxs[l]*multipliers[old2new[l+old_num_row_idxs]];}
	  else {y2+=old_col_idxs[l]*multipliers[old2new[l+old_num_row_idxs]];}
	}
	for (Sparseint p=old.m_array->p[old_c];p<old.m_array->p[old_c+1];++p){ //row pointer
	  //now do the same for the row index, if we haven't decomposed this row before
	  if (!x1x2[old.m_array->i[p]]) { //haven't decomposed this row before?
	    x1x2[old.m_array->i[p]]=new Sparseint[2];
	    std::fill(x1x2[old.m_array->i[p]],x1x2[old.m_array->i[p]]+2,0);
	    index_decompose(old.m_array->i[p],old_num_row_idxs,old_idx_dims.data(),old_row_idxs);
	    for (Sparseint l=0;l<old_num_row_idxs;++l){
	      if (old2new[l]<new_num_row_idxs){x1x2[old.m_array->i[p]][0]+=old_row_idxs[l]*multipliers[old2new[l]];}
	      else {x1x2[old.m_array->i[p]][1]+=old_row_idxs[l]*multipliers[old2new[l]];}
	    }
	  }
	  //check the value is not effectively zero
	  //if (abs(old.m_array->x[p])>SPARSETOL*largest){
	    Sparseint entry;   
	    entry=new_matrix.m_array->nz++;
	    new_matrix.m_array->i[entry]=x1x2[old.m_array->i[p]][1]+y2;
	    new_matrix.m_array->p[entry]=x1x2[old.m_array->i[p]][0]+y1;
	    new_matrix.m_array->x[entry]=conjugate ? conj(old.m_array->x[p]) : old.m_array->x[p];
	    //}
	}
      }
    }
    delete[] old_row_idxs;    
    delete[] old_col_idxs;
    for (Sparseint i=0;i<old.rows();++i){
      delete[] x1x2[i];
    }
    delete[] x1x2;
    delete[] old2new;
    delete[] multipliers;
    return new_matrix.cheap_finalise();
  }
  ////////////////////////////////////////////////////////////////
// a cheaper generalised reshape
  SparseMatrix cheap_reshape(const SparseMatrix& old,Sparseint old_num_row_idxs,Sparseint new_num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<Sparseint>& new_idx_order,const bool conjugate)
  {
    //check finalised
    if (!old.is_finalised()){std::cout << "Input array wasn't finalised?" << std::endl; exit(1);}

    //check product of dims is the same for old
    Sparseint old_prod=1;
    for (std::vector<Sparseint>::const_iterator cit=old_idx_dims.begin();cit!=old_idx_dims.end();++cit){
      old_prod*=(*cit);
    }
    if (old_prod!=old.rows()*old.cols()){std::cout << "Incorrect old params: " << old.rows()*old.cols() << ", " << old_prod << std::endl; exit(1);}
    if (old_idx_dims.size()!=new_idx_order.size()){std::cout << "Incorrect size params, new and old number of indices don't match " << old_idx_dims.size() <<" "<< new_idx_order.size() << std::endl; exit(1);}

    Sparseint old_num_col_idxs=old_idx_dims.size()-old_num_row_idxs;
    Sparseint* multipliers=new Sparseint[old_idx_dims.size()];
    Sparseint new_num_rows=1;
    Sparseint new_num_cols=1;
    Sparseint** x1x2=0;
    Sparseint* old2new=0;
    {
      for (Sparseint i=0;i<new_num_row_idxs;++i){
	multipliers[i]=new_num_rows;
	new_num_rows*=old_idx_dims[new_idx_order[i]];
      }
      for (size_t i=new_num_row_idxs;i<new_idx_order.size();++i){
	multipliers[i]=new_num_cols;
	new_num_cols*=old_idx_dims[new_idx_order[i]];
      }
      {
	x1x2=new Sparseint*[old.rows()];
	std::fill(x1x2,x1x2+old.rows(),static_cast<Sparseint*>(0));
      }
      //make the index mulitpliers
      {
	old2new=new Sparseint[old_idx_dims.size()];
	for (size_t l=0;l<old_idx_dims.size();++l){
	  old2new[new_idx_order[l]]=l;
	}
      }
    }

    SparseMatrix new_matrix(new_num_cols,new_num_rows,old.nz()); //SparseMatrix usually assumes tranposed input...
    //create an array to store the subindices in
    Sparseint* old_col_idxs;
    Sparseint* old_row_idxs;

    old_col_idxs=new Sparseint[old_num_col_idxs];
    old_row_idxs=new Sparseint[old_num_row_idxs];
    for (Sparseint old_c=0;old_c<old.m_array->n;++old_c){ //column
      //take apart the column index if, and only if, the column is not empty
      if(old.m_array->p[old_c]!=old.m_array->p[old_c+1]){
	index_decompose(old_c,old_num_col_idxs,old_idx_dims.data()+old_num_row_idxs,old_col_idxs);
	Sparseint y1=0;
	Sparseint y2=0;
	//and now update old_idxs by their multipliers
	for (Sparseint l=0;l<old_num_col_idxs;++l){
	  if (old2new[l+old_num_row_idxs]<new_num_row_idxs){y1+=old_col_idxs[l]*multipliers[old2new[l+old_num_row_idxs]];}
	  else {y2+=old_col_idxs[l]*multipliers[old2new[l+old_num_row_idxs]];}
	}
	for (Sparseint p=old.m_array->p[old_c];p<old.m_array->p[old_c+1];++p){ //row pointer
	  //now do the same for the row index, if we haven't decomposed this row before
	  if (!x1x2[old.m_array->i[p]]) { //haven't decomposed this row before?
	    x1x2[old.m_array->i[p]]=new Sparseint[2];
	    std::fill(x1x2[old.m_array->i[p]],x1x2[old.m_array->i[p]]+2,0);
	    index_decompose(old.m_array->i[p],old_num_row_idxs,old_idx_dims.data(),old_row_idxs);
	    for (Sparseint l=0;l<old_num_row_idxs;++l){
	      if (old2new[l]<new_num_row_idxs){x1x2[old.m_array->i[p]][0]+=old_row_idxs[l]*multipliers[old2new[l]];}
	      else {x1x2[old.m_array->i[p]][1]+=old_row_idxs[l]*multipliers[old2new[l]];}
	    }
	  }
	  //check the value is not effectively zero
	  Sparseint entry;   
	  entry=new_matrix.m_array->nz++;
	  //doesn't use transpose
	  new_matrix.m_array->i[entry]=x1x2[old.m_array->i[p]][0]+y1;
	  new_matrix.m_array->p[entry]=x1x2[old.m_array->i[p]][1]+y2;
	  new_matrix.m_array->x[entry]=conjugate ? conj(old.m_array->x[p]) : old.m_array->x[p];	  
	}
      }
    }
    delete[] old_row_idxs;    
    delete[] old_col_idxs;
    for (Sparseint i=0;i<old.rows();++i){
      delete[] x1x2[i];
    }
    delete[] x1x2;
    delete[] old2new;
    delete[] multipliers;

    return new_matrix.cheap_no_transpose_finalise();
  }
  //////////////////////////////////////////////////////////
  SparseMatrix permute(const SparseMatrix& A, const std::vector<Sparseint>& column_permutations){
    const Sparseint* q=&column_permutations[0];
    const Sparseint* p=0;
    return SparseMatrix(cs_cl_permute(A.m_array,p,q,1),1);
  }

  /*SparseMatrix permute(const SparseMatrix& A,  const std::vector<Sparseint>& row_permutations, const std::vector<Sparseint>& column_permutations){
    const Sparseint* q=&column_permutations[0];
    const Sparseint* p=&row_permutations[0];
    return SparseMatrix(cs_cl_permute(A.m_array,p,q,1),1);
    }*/

  std::vector<Sparseint> GetNonZeroColumns(const SparseMatrix& M){
    if (!M.is_finalised()){std::cout << "Input array wasn't finalised?" << std::endl; exit(1);}
    std::vector<Sparseint> ans;
    for (Sparseint c=0;c<M.cols();++c){
      if(M.get_p(c)!=M.get_p(c+1)){ans.push_back(c);}
    }
    return ans;
  }

  double SparseMatrix::norm(const std::vector<Sparseint>& cols) const {
    double maxnorm=0.0;
    for (std::vector<Sparseint>::const_iterator col_it=cols.begin();col_it!=cols.end();++col_it){
      double s=0.0;
      for (Sparseint p=this->get_p(*col_it);p<this->get_p(*col_it+1);++p){
	s+=abs(this->get_x(p));
      }
      maxnorm = maxnorm>s ? maxnorm : s;
    }
    return maxnorm;
  }

  double SparseMatrix::sum_column_square_norms(const std::vector<Sparseint>& cols) const {
    double sum=0.0;
    for (std::vector<Sparseint>::const_iterator col_it=cols.begin();col_it!=cols.end();++col_it){
      for (Sparseint p=this->get_p(*col_it);p<this->get_p(*col_it+1);++p){
	sum+=real(conj(this->get_x(p))*(this->get_x(p)));
      }
    }
    return sum;
  }

  double SparseMatrix::square_norm() const {
    double sum=0.0;
    for (Sparseint p=0;p<this->nz();++p){
      sum+=real(conj(this->get_x(p))*(this->get_x(p)));
    }
    return sum;
  }

  SparseSVD SparseMatrix::SVD(const std::vector<std::vector<Sparseint> >& B,size_t D, double min_s_val) const {
    //loop through groups doing svd
#ifndef NDEBUG
    std::cout << "Forming objects for Sparse SVD" << std::endl;
#endif
    std::vector<std::pair<Sparseint,double> > UnsortedValues;
    //U, Vdagger, Values
    //we work with columns of U and columns of Vstar
    SparseMatrix UnsortedU(this->rows(),this->rows());
    SparseMatrix UnsortedVstar(this->cols(),this->cols());
    UnsortedValues.reserve(this->rows() >this->cols() ? this->cols() : this->rows());
    double total_weight(0.0);
    double tolerance=min_s_val < 0.0 ? 0.0 : min_s_val*min_s_val;

#ifndef NDEBUG
    std::cout << "Looping over " << B.size() << " blocks" << std::endl;
#endif
    for (std::vector<std::vector<Sparseint> >::const_iterator cit=B.begin();cit!=B.end();++cit){
#ifndef NDEBUG
      std::cout << "Block " <<  std::distance(B.begin(),cit+1) << std::endl;
      //check block isn't effectively empty?
      std::cout << "Calculating block weight" << std::endl;
#endif
      double block_weight(sum_column_square_norms(*cit));
      total_weight+=block_weight;
      if (block_weight<= tolerance){
#ifndef NDEBUG
	std::cout << "This block is empty, skipping" << std::endl; 
#endif
	continue;
      } //skip block
      std::vector<Sparseint> cols;
#ifndef NDEBUG
      std::cout << "Getting non empty block columns" << std::endl;
#endif
      for (std::vector<Sparseint>::const_iterator col_it=cit->begin();col_it!=cit->end();++col_it){
  	if (this->get_p(*col_it)!=this->get_p(*col_it+1)){
  	  cols.push_back(*col_it);
  	}
      }
#ifndef NDEBUG
      std::cout << "Forming Dense Translation block" << std::endl;
#endif
      TranslationBlock<DenseMatrix> TB(*this,cols /* *cit*/ ); //form block
#ifndef NDEBUG
      std::cout << "Performing Dense SVD" << std::endl;
#endif
      DenseSVD Blockans(std::move(TB.Block.SVD())); //do svd on block (which will destroy the original block)
      // now read through block answer and feed back into sparse
#ifndef NDEBUG
      std::cout << "Reading back into Sparse" << std::endl;
#endif
      for (Sparseint v=0;v<Blockans.ValuesSize();++v){
	//read column of U into sparse
	for (Sparseint i=0;i<Blockans.U.rows();++i){
	  //because cols of U are normalised vectors, and these vectors are not absurdly long, we drop the smallest values.
	  if (abs(Blockans.U.Value(i,v))>SPARSETOL){
	    UnsortedU.entry(TB.RowLookups[i],UnsortedValues.size(),Blockans.U.Value(i,v));
	  }
	}
	//read row of Vdagger into sparse 
	for (Sparseint j=0;j<Blockans.Vdagger.cols();++j){
	  if (abs(Blockans.Vdagger.Value(v,j))>SPARSETOL){
	    UnsortedVstar.entry(TB.ColLookups[j],UnsortedValues.size(),Blockans.Vdagger.Value(v,j));
	  }
	}
	UnsortedValues.push_back(std::pair<Sparseint,double>(UnsortedValues.size(),Blockans.Values[v]));
      }
    }

    if (UnsortedValues.size()==0){
      std::cout << "Error: No singular values. Empty blocks?" << std::endl;
      this->print();
      for (auto&& bc : B){
	std::cout << "Block column indices: ";
	for (auto c : bc){
	  std::cout << c << " ";
	}
	std::cout << std::endl;
      }
      exit(1);
    }

#ifndef NDEBUG
    std::cout << "Sorting " << UnsortedValues.size()  << " singular values" << std::endl;
#endif
    //we will always sort, in case we want to truncate later
    std::sort(UnsortedValues.begin(),UnsortedValues.end(),singular_value_compare);
    size_t length=(D>0 && D< UnsortedValues.size()) ? D : UnsortedValues.size();

    //We should also probably reject singular vals that are smaller the M_EPS times the largest singular value...
    std::vector<Sparseint> sortindices;
    sortindices.reserve(length);
    sortindices.push_back(UnsortedValues[0].first);
    std::vector<double> Values;
    Values.reserve(length);
    Values.push_back(UnsortedValues[0].second);
    double kept_weight(UnsortedValues[0].second*UnsortedValues[0].second);
    for (size_t s=1;s<length;++s){
      if (UnsortedValues[s].second < SPARSETOL*UnsortedValues.begin()->second) {
	length=s;
	std::cout << "Truncating further due to very small singular values." <<std::endl; 
	break;
      }
      sortindices.push_back(UnsortedValues[s].first);
      kept_weight+=UnsortedValues[s].second*UnsortedValues[s].second;
      Values.push_back(UnsortedValues[s].second);
    }

    double discarded_weight(0.0);
    for (size_t s=length;s<UnsortedValues.size();++s){
      discarded_weight+=UnsortedValues[s].second*UnsortedValues[s].second;
    }

#ifndef NDEBUG
    std::cout << "SVD()" << std::endl;
    std::cout << "Largest s val: " << UnsortedValues[0].second << ", Smallest s val: " << UnsortedValues[length-1].second << std::endl <<std::endl;
    std::cout << "Total weight is: " << total_weight << std::endl;
    std::cout << "Kept weight is: " << kept_weight << std::endl;
    std::cout << "Discarded weight is: " << discarded_weight << std::endl <<std::endl;
    std::cout << "Resizing arrays" << std::endl;
#endif
    //Unsorted are now in fact sorted
    UnsortedU.resize_finalise(this->rows(),UnsortedValues.size()); //get rid of empty cols
    UnsortedVstar.resize_finalise(this->cols(),UnsortedValues.size());
#ifndef NDEBUG
    std::cout << "Populating output object" << std::endl;
#endif
    return SparseSVD(std::move(Values),std::move(UnsortedU.ExtractColumns(sortindices)),std::move(UnsortedVstar.ExtractColumns(sortindices).transpose()),kept_weight,discarded_weight);
  }
  //////////////////////////////////////////////////////////
  std::vector<double> SparseMatrix::SVD() const{
    //make a DenseMatrix
    return (this->to_dense()).SingularValues();
  }
  //////////////////////////////////////////////////////////
  //first variant finds all eigenvalues and vectors using blocks and a dense method
  SparseHED SparseMatrix::HED(const std::vector<std::vector<Sparseint> >& B) const{
    //std::cout << "METHOD 1" << std::endl;
    if (!this->is_finalised()){std::cout << "Not finalised!" << std::endl; exit(1);}
    //std::cout << "HED on matrix, size " << this->rows() << " * " << this->cols() << std::endl;
    if (this->rows()!=this->cols()){std::cout << "Hermitian matrix not square!" << std::endl; exit(1);}
    SparseHED ans(this->rows(),this->cols());
    for (std::vector<std::vector<Sparseint> >::const_iterator cit=B.begin();cit!=B.end();++cit){
      TranslationBlock<DenseMatrix> TB(*this,*cit,*cit); //form block
      DenseHED Blockans=TB.Block.HED(); //do eigen decomposition on block (which will destroy the original block)
      // now read through block answer and feed back into sparse
      //copy across eigenvalues
      //copy in eigenvectors
      for (Sparseint v=0;v<Blockans.ValuesSize();++v){
	//read column of EigenVectors into sparse
	for (Sparseint i=0;i<Blockans.EigenVectors.rows();++i){
	  if (abs(Blockans.EigenVectors.Value(i,v))>SPARSETOL){
	    //ans.EigenVectors.entry(TB.ColLookups[i],ans.ValuesSize(),Blockans.EigenVectors.Value(i,v));
	    ans.EigenVectors.entry(TB.ColLookups[i],ans.ValuesSize(),Blockans.EigenVectors.Value(i,v));
	  }
	}
	ans.Values.push_back(Blockans.Values[v]);
      }   
    }
    //at this point the sparsematrices are still in triplet form
    //so need to drop any zero rows (from Vdagger) or cols from (U).
    ans.EigenVectors.finalise();//resize_finalise(this->rows(),ans.Values.size());
    return ans;
  }
  //////////////////////////////////////////////////////////
  //second variant finds some eigenvalues and vectors for a block using dense or sparse methods
  SparseHED SparseMatrix::HED(const std::vector<Sparseint>& B, Sparseint numevals, char which[3],SparseMatrix* initial) const{
    //std::cout << "METHOD 2: "<<std::endl;
    if (!this->is_finalised()){std::cout << "Not finalised!" << std::endl; exit(1);}
    if (this->rows()!=this->cols()){std::cout << "Hermitian matrix not square!" << std::endl; exit(1);}
    //if (numevals>=this->rows()){std::cout << "Requesting too many eigenvalues and vectors " << numevals << ", " << this->rows() << std::endl; exit(1);}
    //now choose dense or sparse method for few eigenvals
    //if we have matrices smaller than 400*400 then dense is probably fine
    //larger matrices should use a sparse method
    if (this->rows()<=400){ //use dense storage for block
      //std::cout << "Method 2 with dense storage" <<std::endl;
      TranslationBlock<DenseMatrix> TB(*this,B,B); //form block
      
      numevals=numevals>TB.Block.rows() ? TB.Block.rows() : numevals;
      SparseHED ans(this->rows(),numevals);

      DenseHED Blockans=TB.Block.HED(numevals,which);
      // now read through block answer and feed back into sparse
      //copy across eigenvalues
      //copy in eigenvectors
      for (Sparseint v=0;v<Blockans.ValuesSize();++v){
	//read column of EigenVectors into sparse
	for (Sparseint i=0;i<Blockans.EigenVectors.rows();++i){
	  if (abs(Blockans.EigenVectors.Value(i,v))>SPARSETOL){
	      //ans.EigenVectors.entry(TB.ColLookups[i],ans.ValuesSize(),Blockans.EigenVectors.Value(i,v));
	    ans.EigenVectors.entry(TB.ColLookups[i],ans.ValuesSize(),Blockans.EigenVectors.Value(i,v));
	  }
	}
	ans.Values.push_back(Blockans.Values[v]);
      }

      ans.EigenVectors.finalise();
      return ans;

    }
    else { //use sparse storage for block
      TranslationBlock<SparseMatrix> SpTB(*this,B,B); //form block
      numevals=numevals>SpTB.Block.rows() ? SpTB.Block.rows() : numevals;
      SparseHED ans(this->rows(),numevals);

      SparseMatrix* Guessptr;
      SparseMatrix Guess;
      if (initial) {Guess=SpTB.ReverseTranslateRows(*initial); Guessptr=&Guess;std::cout << "Allocating Guess" <<std::endl;}
      else {Guessptr=nullptr;}
      //now calls eigensolver on just the translation block.
      //if translation block is smaller than 100 x 100, it will end up using a dense method.
      SparseHED SpBlockans(SpTB.Block.HED(numevals,which,Guessptr));
      ans.Values=SpBlockans.Values; //copy eigenvalues over
      for (Sparseint c=0;c<numevals;++c){
	//ans.Values.push_back(SpBlockans.Values.at(c));
	for (Sparseint p=SpBlockans.EigenVectors.get_p(c);p<SpBlockans.EigenVectors.get_p(c+1);++p){
	  ans.EigenVectors.entry(SpTB.ColLookups[SpBlockans.EigenVectors.get_i(p)],c,SpBlockans.EigenVectors.get_x(p));
	}
      }
      
      ans.EigenVectors.finalise();
      return ans;
    }
  }

  //////////////////////////////////////////////////////////
  //third variant, find some eigenvals and vectors, without further blocking via dense method or arpack
  //used by the second variant above
  //includes some massaging because arpack doesn't assume Hermitian form
  SparseHED SparseMatrix::HED(Sparseint numevals, char which[3],SparseMatrix* initial) const {
    //std::cout << "METHOD 3" << std::endl;
    if (this->rows()!=this->cols()){std::cout << "Matrix not square for HED!" << std::endl; exit(1);}
    numevals=numevals>this->rows() ? this->rows() : numevals;
    SparseHED ans(this->rows(),numevals);
    std::complex<double>* Evecs = new std::complex<double>[this->rows()*numevals];
    std::complex<double>* Evals = new std::complex<double>[numevals];
    //can we guarantee arpack won't mess the array up?
    bool arpack_error=1;
    if (this->rows()>100){
      arpack_error=arpack::cpparpack(this->m_array,this->rows(),numevals,Evals, Evecs, which,initial ? initial->m_array : NULL);
    }
    if (arpack_error){ 
      //use dense storage for block
      //if arpack returns an error we have to use a dense method
      //std::cout << "Method 3 using dense storage..." <<std::endl;
     DenseHED densedecomp(this->to_dense().HED(numevals,which));
      std::copy(densedecomp.Values,densedecomp.Values+numevals,Evals);
      for (Sparseint j=0;j<numevals;++j){
	for (Sparseint i=0;i<this->rows();++i){
	  Evecs[i+j*this->rows()]=densedecomp.EigenVectors.Value(i,j);
	}
      }
    }
    for (size_t v=0;v<static_cast<size_t>(numevals);++v){
      std::cout << "Eval: " << Evals[v] <<std::endl;
      ans.Values.push_back(Evals[v].real()); //arpack doesn't assume a Hermitian form
      double maxabsvalue=0.0;
      std::complex<double> maxvalue=0.0;
      //first pass find largest abs value
      for (Sparseint i=0;i<this->rows();++i){
	if (abs(Evecs[i+v*this->rows()])>maxabsvalue){
	  maxabsvalue = abs(maxvalue=Evecs[i+v*this->rows()]);
	}
      }
      maxvalue=maxabsvalue/maxvalue;
      //now get rid of possible phase
      for(Sparseint i=0;i<this->rows();++i){
	std::complex<double> value=Evecs[i+v*this->rows()]*maxvalue;
	if (abs(value)>SPARSETOL*maxabsvalue){
	  ans.EigenVectors.entry(i,v,FixComplexPrecision(value));
	}
      }
    }
    ans.EigenVectors.finalise();
    delete[] Evecs;
    delete[] Evals;
    return ans;
  }

  SparseED SparseMatrix::RightED(const std::vector<Sparseint>& RowB,const std::vector<Sparseint>& ColB, Sparseint numevals, char which[3],SparseMatrix* initial) const{
    if (!this->is_finalised()){std::cout << "Not finalised!" << std::endl; exit(1);}
    if (this->rows()!=this->cols()){std::cout << "Matrix not square for RightED!" << std::endl; exit(1);}
    SparseED ans(this->cols(),numevals);
    //use sparse storage for block
    TranslationBlock<SparseMatrix> SpTB(*this,RowB,ColB); //form block
    SparseMatrix* Guessptr;
    SparseMatrix Guess;
    if (initial) {Guess=SpTB.ReverseTranslateRows(*initial); Guessptr=&Guess;}
    else {Guessptr=NULL;}
    SparseED SpBlockans(SpTB.Block.ED(numevals,which,Guessptr));
    ans.Values=SpBlockans.Values; //copy eigenvalues over
    for (Sparseint c=0;c<numevals;++c){
      //ans.Values.push_back(SpBlockans.Values.at(c));
      for (Sparseint p=SpBlockans.EigenVectors.get_p(c);p<SpBlockans.EigenVectors.get_p(c+1);++p){
	ans.EigenVectors.entry(SpTB.ColLookups[SpBlockans.EigenVectors.get_i(p)],c,SpBlockans.EigenVectors.get_x(p));
      }
    }
    ans.EigenVectors.finalise();
    return ans;
  }

  SparseED SparseMatrix::RightED(Sparseint numevals, char which[3],SparseMatrix* initial) const{
    if (!this->is_finalised()){std::cout << "Not finalised!" << std::endl; exit(1);}
    if (this->rows()!=this->cols()){std::cout << "Matrix not square for RightED!" << std::endl; exit(1);}
    //output will be left eigenvectors, stored as column vectors
    return this->ED(numevals,which,initial);
  }

  //getting the left eigenvectors is a bit trickier
  //requires a transpose of the translation block, use of arpack then transpose back
  SparseED SparseMatrix::LeftED(const std::vector<Sparseint>& RowB,const std::vector<Sparseint>& ColB, Sparseint numevals, char which[3],SparseMatrix* initial) const{
    if (!this->is_finalised()){std::cout << "Not finalised!" << std::endl; exit(1);}
    if (this->rows()!=this->cols()){std::cout << "Matrix not square for LeftED!" << std::endl; exit(1);}
    //output will be left eigenvectors, stored as column vectors
    SparseED ans(this->rows(),numevals);
    //use sparse storage for block
    TranslationBlock<SparseMatrix> SpTB(*this,RowB,ColB); //form block
    //but we need the transpose
    SparseMatrix* Guessptr;
    SparseMatrix Guess;
    //initial guess (if used) should still be fed in as a column vector
    if (initial) {Guess=SpTB.ReverseTranslateRows(*initial); Guessptr=&Guess;}
    else {Guessptr=NULL;}
    //below does a transpose, messing up the translation block's sparsematrix
    SparseED SpBlockans(SpTB.Block.transpose().ED(numevals,which,Guessptr));
    ans.Values=SpBlockans.Values; //copy eigenvalues over
    for (Sparseint c=0;c<numevals;++c){
      for (Sparseint p=SpBlockans.EigenVectors.get_p(c);p<SpBlockans.EigenVectors.get_p(c+1);++p){
	ans.EigenVectors.entry(SpTB.RowLookups[SpBlockans.EigenVectors.get_i(p)],c,SpBlockans.EigenVectors.get_x(p));
      }
    }
    ans.EigenVectors.finalise();
    return ans;
  }

  SparseED SparseMatrix::LeftED(Sparseint numevals, char which[3],SparseMatrix* initial) const{
    if (!this->is_finalised()){std::cout << "Not finalised!" << std::endl; exit(1);}
    if (this->rows()!=this->cols()){std::cout << "Matrix not square for LeftED!" << std::endl; exit(1);}
    //output will be left eigenvectors, stored as column vectors
    return this->copy_transpose().ED(numevals,which,initial);
  }

  SparseED SparseMatrix::ED(Sparseint requested_numevals, char which[3],SparseMatrix* initial) const {
    if (this->rows()!=this->cols()){std::cout << "Matrix not square for ED! " << this->rows() << " " << this->cols() << std::endl; exit(1);}
    std::complex<double>* Evecs = new std::complex<double>[this->rows()*requested_numevals];
    std::complex<double>* Evals = new std::complex<double>[requested_numevals];
    SparseED ans(this->rows(),requested_numevals);

    //Run arpack and check for errors
    if(this->cols() < 10 || arpack::cpparpack(this->m_array,this->rows(),requested_numevals,Evals, Evecs, which,initial ? initial->m_array : NULL)){
      //std::cout << "Arpack error for array size " << this->rows() << " " << this->cols() << std::endl;
      if (this->rows()==1 && this->cols()==1){
	Evecs[0]=1.0; Evals[0]=this->m_array->x[0];
	std::cout << "Trivial answer" << std::endl;
      }
      else {
//turn initial to dense
	std::cout << "Using dense method for Non Hermitian eigenvalues..." << std::endl;

	std::complex<double>* dEvecs= new std::complex<double>[this->rows()*this->rows()];
	std::complex<double>* dEvals = new std::complex<double>[this->rows()];

	std::complex<double>* dm=new std::complex<double>[this->rows()*this->rows()];
	this->to_dense(dm);
	densefuncs::diagonalise_with_lapack_nh(this->rows(), dm, dEvecs, dEvals);
	delete[] dm;

	//now sort
	double weight=0.0;
	Sparseint largest=0;
	for (Sparseint i=0; i<this->rows();++i){
	  if (abs(dEvals[i])>weight) {largest=i;weight=abs(dEvals[i]);}
	}
	if (largest!=0){
	  std::swap_ranges(dEvals,dEvals+1,&dEvals[largest]);
	  std::swap_ranges(dEvecs,dEvecs+this->rows(),&(dEvecs[largest]));
	}
	//copy into
	std::copy(dEvals,dEvals+requested_numevals,Evals);
	std::copy(dEvecs,dEvecs+this->rows()*requested_numevals,Evecs);
	delete[] dEvals;
	delete[] dEvecs;
      }
    }
    for (size_t v=0;v<static_cast<size_t>(requested_numevals);++v){
      ans.Values.push_back(FixComplexPrecision(Evals[v]));
      double maxabsvalue=0.0;
      std::complex<double> maxvalue=0.0;
      //first pass find largest abs value
      for (Sparseint i=0;i<this->rows();++i){
	if (abs(Evecs[i+v*this->rows()])>maxabsvalue){
	  maxabsvalue = abs(maxvalue=Evecs[i+v*this->rows()]);
	}
      }
      maxvalue=maxabsvalue/maxvalue;
      //now get rid of possible phase
      for(Sparseint i=0;i<this->rows();++i){
	std::complex<double> value=Evecs[i+v*this->rows()]*maxvalue;
	if (abs(value)>SPARSETOL*maxabsvalue){
	  ans.EigenVectors.entry(i,v,FixComplexPrecision(value));
	}
      }
    }
    ans.EigenVectors.finalise();
    delete[] Evecs;
    delete[] Evals;
    return ans;
  }

  std::pair<SparseMatrix,SparseMatrix> SparseHermitianDecomposition(DenseHED decomp){
    std::pair<SparseMatrix,SparseMatrix> ans(SparseMatrix(decomp.EigenVectors.cols(),decomp.EigenVectors.rows()),SparseMatrix(decomp.EigenVectors.rows(),decomp.EigenVectors.cols()));
      //SparseMatrix ans(decomp.EigenVectors.cols(),decomp.EigenVectors.rows());
      //SparseMatrix ansinv(decomp.EigenVectors.rows(),decomp.EigenVectors.cols());
    for (Sparseint j=0;j<decomp.EigenVectors.cols();++j){
      for (Sparseint i=0;i<decomp.EigenVectors.rows();++i){
	if (decomp.Values[j]<=0.0){std::cout << "Matrix isn't positive, can't do Hermitian decomp" << std::endl; exit(1);}
	{
	  //complex<double> element=sqrt(decomp.Values[j])*decomp.EigenVectors.Value(i,j);
	  if (abs(decomp.EigenVectors.Value(i,j))>SPARSETOL){
	    ans.first.entry(j,i,sqrt(decomp.Values[j])*decomp.EigenVectors.Value(i,j));
	    ans.second.entry(i,j,decomp.EigenVectors.Value(i,j)/sqrt(decomp.Values[j]));
	  }
	}
      }
    }
    ans.first.finalise();
    ans.second.finalise();
    return ans;
  }

  SparseMatrix Exponentiate(const SparseHED& decomp,std::complex<double> factor){
    std::vector<std::complex<double> > diagpart;
    for (std::vector<double>::const_iterator cit=decomp.Values.begin();cit!=decomp.Values.end();++cit ){
      diagpart.push_back(exp(factor*(*cit)));
    }
    return (decomp.EigenVectors*SparseMatrix(diagpart))*decomp.EigenVectors.copy_dagger();
  }

  void DumbExtractWithZeros(const SparseMatrix& S,const std::vector<Sparseint>& Allowed, std::complex<double>* out){//assumes the entries in allowed are ordered
    if (S.cols()!=1) {std::cout << "Not a sparse VECTOR!" << std::endl; exit(1);}
    uMPXInt p=0, a=0;
    //std::cout << "Svec length " << S.rows() << " Allowed length " << Allowed.size() << std::endl;
    for (std::vector<Sparseint>::const_iterator cit=Allowed.begin();cit!=Allowed.end();++cit){//increment 'Allowed' and 'a'
      if (S.get_i(p)==*cit){//does sparse row == next allowed entry row?
	out[a++]=S.get_x(p++); //if yes, populated answer, move to next non zero
      }
      else {out[a++]=0.0;} //else pad with zero
      if (p>=S.get_p(1)){//if we've gone through all the non zeros, then finish off with zeros and exit loop
	while (a<Allowed.size()){
	  out[a++]=0.0;
	}
	break;
      }
    }
  }

  void DumbExtractUpdate(const SparseMatrix& S,const std::vector<Sparseint>& Allowed, std::complex<double>* out, std::complex<double> weight){//assumes the entries in allowed are ordered
    if (S.cols()!=1) {std::cout << "Not a sparse VECTOR!" << std::endl; exit(1);}
    uMPXInt p=0, a=0;
    for (std::vector<Sparseint>::const_iterator cit=Allowed.begin();cit!=Allowed.end();++cit,++a){//increment 'Allowed' and 'a'
      if (S.get_i(p)==*cit){//does sparse row == next allowed entry row?
	out[a]+=weight*S.get_x(p++); //if yes, populated answer, move to next non zero
      }
      
      if (p>=S.get_p(1)){//if we've gone through all the non zeros, then exit loop
	break;
      }
    }
  }

  DenseMatrix ExtractWithZeros(const SparseMatrix& S,const std::vector<Sparseint>& Allowed){
    if (S.cols()!=1) {std::cout << "Not a sparse VECTOR!" << std::endl; exit(1);}
    std::complex<double>* array=new std::complex<double>[S.rows()];
    DumbExtractWithZeros(S,Allowed,array); //fill array
    return DenseMatrix(S.rows(),1,&array);
  }

  void order_rows(cs_cl* U){
    if (!U){std::cout << "Error" <<std::endl; exit(1);}
    if (U->nz!=-1){std::cout << "Error" <<std::endl; exit(1);}
    if (U->m!=1){
    //record original col counts
    SuiteSparse_long* Cp = U->p;
    SuiteSparse_long* Ci = U->i;
    std::complex<double>* Cx = U->x ;
    cs_cl* A=cs_cl_transpose(U,-1);
    //order by double transpose, without recounting rows
    SuiteSparse_long m = A->m; 
    SuiteSparse_long n = A->n; 
    SuiteSparse_long* Ap = A->p;
    SuiteSparse_long* Ai = A->i; 
    std::complex<double>* Ax = A->x;
    SuiteSparse_long* w=(SuiteSparse_long*)cs_malloc(m,sizeof(SuiteSparse_long));
    std::copy(Cp,Cp+m,w);//need a copy in w
    for (SuiteSparse_long j = 0; j < n; j++){
      for (SuiteSparse_long p = Ap[j]; p < Ap[j+1]; p++){
	SuiteSparse_long q;
	Ci[q = w [Ai[p]]++] = j; /* place ATranspose(i,j) as entry A(j,i) */
	Cx[q] = Ax[p];
      }
    }
    cs_cl_spfree(A);
    cs_free(w);
    }
    //return A;
  }

  void order_rows_and_feed_to_answer(cs_cl* U, cs_cl* P,const SuiteSparse_long Pcol){//destroys unsorted U as it goes along, Pcol is the insertion column
    //record original col counts
    cs_cl* UTranspose=cs_cl_transpose(U,-1); //orders, but gives us the transpose
    //order by double transpose, without recounting rows
    SuiteSparse_long m = UTranspose->m; 
    SuiteSparse_long n = UTranspose->n; 
    SuiteSparse_long* UTp = UTranspose->p;
    SuiteSparse_long* UTi = UTranspose->i; 
    std::complex<double>* UTx = UTranspose->x;
    //P already allocated..., need to start at last nonzero of P and update
    SuiteSparse_long* Pi = P->i;
    std::complex<double>* Px = P->x;
    //to be thread safe must use a separate counter for nonzeros of P!!!!
    cs_free(U->i);  U->i=nullptr;
    cs_free(U->x);  U->x=nullptr;
    SuiteSparse_long* w = U->p;  U->p=nullptr;
    U->nzmax=0;
    cs_cl_spfree(U); //have hung on to U->p.
    //P->p should have been set already
    SuiteSparse_long offset=P->p[Pcol];
    for (SuiteSparse_long j = 0; j < n; j++){
      for (SuiteSparse_long p = UTp[j]; p < UTp[j+1]; p++){
	SuiteSparse_long q= offset+w[UTi[p]]++;
	//P->p already exists, so doesn't need copying in
	//P->p[Pcol] should give the correct starting position
	Pi[q] = j; //use the offset from P
	Px[q] = UTx[p];
	//w[UTi[p]]++;//increment pointer each time we add a value
      }
    }
    cs_free(w); //finally deallocated old U->p = w
    cs_cl_spfree(UTranspose);
  }

  void feed_to_answer(cs_cl* U, cs_cl* P,const SuiteSparse_long Pcol){//destroys unsorted U as it goes along, Pcol is the insertion column
    SuiteSparse_long Unz=U->p[U->n];
    std::copy(U->i,U->i+Unz,P->i+P->p[Pcol]);
    std::copy(U->x,U->x+Unz,P->x+P->p[Pcol]);
    cs_cl_spfree(U);
  }

  cs_cl* multiply_by_blocked_cols(const cs_cl* A, const cs_cl*B, SuiteSparse_long start, SuiteSparse_long end, const SuiteSparse_long* ans_col_ptrs){
    cs_cl *C=nullptr;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL);      /* check inputs */
    if (A->n != B->m) return (NULL);
    SuiteSparse_long nzmax=ans_col_ptrs[end]-ans_col_ptrs[start];
    SuiteSparse_long nz=0;
    SuiteSparse_long m = A->m;
    //SuiteSparse_long anz = A->p [A->n];
    //SuiteSparse_long n = B->n; 
    SuiteSparse_long* Bp = B->p;
    SuiteSparse_long* Bi = B->i ;
    std::complex<double>* Bx = B->x;
    //SuiteSparse_long bnz = Bp [B->n];
    SuiteSparse_long* w = (SuiteSparse_long*)cs_calloc (m, sizeof (SuiteSparse_long));
    std::complex<double>* x = (std::complex<double>*)cs_malloc (m, sizeof (std::complex<double>));
    C = cs_cl_spalloc (m, end-start, nzmax, 1, 0); //we have already figured out the number of nonzeros, so can do the correct allocation now       
    SuiteSparse_long* Cp = C->p;
    for (SuiteSparse_long j = start; j < end; j++){
      SuiteSparse_long* Ci = C->i; 
      std::complex<double>* Cx = C->x;
      Cp[j-start] = nz;                  /* column j of C starts here */
      for (SuiteSparse_long p = Bp[j]; p < Bp[j+1]; p++){
	nz = cs_cl_scatter (A, Bi[p], Bx[p], w, x, j+1, C, nz);
      }
      for (SuiteSparse_long p = Cp[j-start]; p < nz; p++) Cx[p] = x [Ci[p]];
    }
    Cp [end-start] = nz ;                       /* finalize the last column of C */
    assert(nz==nzmax);
    cs_free(w);
    cs_free(x);
    return C;    
  }

#if defined(USETBB)
  struct reduce_nz {
    SuiteSparse_long nz_; // accumulating integer
    const cs_cl* lhs_;
    const cs_cl* rhs_;
    SuiteSparse_long* ans_col_ptrs_;
    reduce_nz(const cs_cl* lhs,const cs_cl* rhs)  : lhs_(lhs),rhs_(rhs),ans_col_ptrs_(nullptr) {
      nz_=0; //init
      ans_col_ptrs_=(SuiteSparse_long*) cs_malloc(rhs_->n+1,sizeof(SuiteSparse_long)); 
      ans_col_ptrs_[0]=0;
    }
    // splitting constructor required by TBB
    reduce_nz( reduce_nz& rb, tbb::split ) : lhs_(rb.lhs_), rhs_(rb.rhs_),ans_col_ptrs_(rb.ans_col_ptrs_){
      nz_=0;
    }
    // the main computation method
    void operator()(const tbb::blocked_range<SuiteSparse_long>& r) { 
      // closely resembles the original serial loop
      int* Flag = new int[lhs_->m];
      std::fill(Flag,Flag+lhs_->m,0);
      std::vector<SuiteSparse_long> non_zero_rows;
      for (SuiteSparse_long col=r.begin(); col<r.end(); ++col){ // iterates over a subrange
	SuiteSparse_long nz_in_this_col=0;
	for (SuiteSparse_long rp=rhs_->p[col];rp<rhs_->p[col+1];++rp){
	  SuiteSparse_long j=rhs_->i[rp];
	  for (SuiteSparse_long lp=lhs_->p[j];lp<lhs_->p[j+1];++lp){
	    SuiteSparse_long i=lhs_->i[lp];
	    if (!Flag[i]){
	      Flag[i]=1;
	      ++nz_;
	      ++nz_in_this_col;
	      non_zero_rows.push_back(i); //record these to reset at end of cycle
	    }
	  }
	}
	ans_col_ptrs_[col+1]=nz_in_this_col;//thread safe because every col writes to a different memory location
	//only reset the values we need to
	for (auto element : non_zero_rows){
	  Flag[element]=0;
	}
	non_zero_rows.clear(); //faster push back next time through, because less allocations needs to occur.
      }
      delete[] Flag;
    }
    // the method to reduce computations accumulated in two bodies,
    // unnecessary really, but provides a check.
    void join( reduce_nz& rb ) {
      nz_+=rb.nz_;
    }
  };

  struct reduce_sparse {
    const cs_cl* lhs_;
    const cs_cl* rhs_;
    bool NoSort_;
    cs_cl* ans_;
    reduce_sparse(const cs_cl* lhs,const cs_cl* rhs,SuiteSparse_long* ans_col_ptrs,bool NoSort)  : lhs_(lhs),rhs_(rhs),NoSort_(NoSort),ans_(nullptr){
      ans_ = (cs_cl*)cs_cl_calloc(1, sizeof (cs_cl)); 
      ans_->m = lhs_->m;                             
      ans_->n = rhs_->n;
      ans_->nzmax = ans_col_ptrs[rhs_->n];
      ans_->nz = -1;            
      ans_->p = ans_col_ptrs;
      ans_->i = (SuiteSparse_long*)cs_malloc (ans_->nzmax, sizeof (SuiteSparse_long));
      ans_->x = (std::complex<double>*)cs_malloc (ans_->nzmax, sizeof (std::complex<double>));
    }
    reduce_sparse(reduce_sparse& rb, tbb::split) : lhs_(rb.lhs_), rhs_(rb.rhs_),NoSort_(rb.NoSort_),ans_(rb.ans_){
    }
    void operator()(const tbb::blocked_range<SuiteSparse_long>& r) { 
      SuiteSparse_long start_col=r.begin();
      SuiteSparse_long end_col;
      for (SuiteSparse_long q=r.begin();q!=r.end();++q) end_col=q+1; //just in case r.end() doesn't have epxected behaviour

      cs_cl* local_answer(multiply_by_blocked_cols(lhs_,rhs_,start_col,end_col,ans_->p));//steals ans_->p to make new answer
      //feed_to_answer(local_answer,ans_,start_col);//moves to full answer and destroys local answer
      if (NoSort_){
	feed_to_answer(local_answer,ans_,start_col);//moves to full answer and destroys local answer
      }
      else {
	order_rows_and_feed_to_answer(local_answer,ans_,start_col);//sorts rows and destroys local answer
      }
      local_answer=nullptr;
    }
    void join(reduce_sparse& rb){}
  };
#endif

  Sparseint* allocate_nonzero_totals(const cs_cl* A, const cs_cl* B){
    if (B->p[B->n]==0 || A->p[A->n]==0){ //epmty array!!
      SuiteSparse_long* ans_col_ptrs_=(SuiteSparse_long*) cs_malloc(B->n+1,sizeof(SuiteSparse_long)); 
      std::fill(ans_col_ptrs_,ans_col_ptrs_+B->n+1,0.0);
      return ans_col_ptrs_;
    }
    else {
#if defined(USETBB)
      Sparseint guess_chunksize=TBBNZ*B->n/(B->p[B->n]);
      if (guess_chunksize<B->n){
	//std::cout << "THREADED" <<std::endl;
	reduce_nz nz_body(A,B);
	tbb::parallel_reduce(tbb::blocked_range<SuiteSparse_long>(0,B->n), nz_body);
	//std::cout << 10000*B->n/(B->p[B->n]) << std::endl;
	for (Sparseint c=1;c<=B->n;++c){//cumulative sum
	  nz_body.ans_col_ptrs_[c]+=nz_body.ans_col_ptrs_[c-1];
	}
	return nz_body.ans_col_ptrs_;
      }
#endif
      //std::cout << "NOT THREADED" <<std::endl;
      SuiteSparse_long* ans_col_ptrs_=(SuiteSparse_long*) cs_malloc(B->n+1,sizeof(SuiteSparse_long)); 
      ans_col_ptrs_[0]=0;
      int* Flag = new int[A->m];
      std::fill(Flag,Flag+A->m,0);
      std::vector<SuiteSparse_long> non_zero_rows;
      for (SuiteSparse_long col=0; col<B->n; ++col){ // iterates over a subrange
	SuiteSparse_long nz_in_this_col=0;
	for (SuiteSparse_long rp=B->p[col];rp<B->p[col+1];++rp){
	  SuiteSparse_long j=B->i[rp];
	  for (SuiteSparse_long lp=A->p[j];lp<A->p[j+1];++lp){
	    SuiteSparse_long i=A->i[lp];
	    if (!Flag[i]){
	      Flag[i]=1;
	      ++nz_in_this_col;
	      non_zero_rows.push_back(i);
	    }
	  }
	}
	ans_col_ptrs_[col+1]=nz_in_this_col;//thread safe because every col writes to a different memory location
	//only reset the values we need to
	for (auto element : non_zero_rows){
	  Flag[element]=0;
	}
	non_zero_rows.clear(); //faster push back next time through, because less allocations needs to occur.
      }
      delete[] Flag;
      for (Sparseint c=1;c<=B->n;++c){//cumulative sum
	ans_col_ptrs_[c]+=ans_col_ptrs_[c-1];
      }
      return ans_col_ptrs_;
    }
  }

  SparseMatrix sparse_multiply_ajaj(const SparseMatrix& lhs, const SparseMatrix& rhs, bool NoSort){
    //first figure out the number of nonzeros in output
    //this is cleaner in terms of reallocations, but might be slower.
    //it should be easy to thread though as we are just accumulating (parallel reducing) nz
    //need an array for flagging whether an entry appears in ans or not.
    if (lhs.cols()!=rhs.rows()) {
      std::cout << "Invalid dimensions for matrix multiplication"<< std::endl;
      std::cout << lhs.cols() << " " << rhs.rows() <<std::endl;
      //exit(1);
    }
    //Threading is over cols of rhs, so that needs to be a decent size.
    //However if number of non zeros is really small, then that is probably not worth it either

    //Are both matrices dense? If so, use dense method.
    if (lhs.is_dense() && rhs.is_dense()){
      return dense_dense_multiply(lhs,rhs);
    }

    //first calculate the number of nonzeros in the answer
    Sparseint* col_ptrs=allocate_nonzero_totals(lhs.m_array,rhs.m_array);// col_ptrs[rhs.cols()]=number of non zeros in answer
    if (col_ptrs[rhs.cols()]==0){//empty
      cs_free(col_ptrs);
      return std::move(SparseMatrix(lhs.rows(),rhs.cols(),1).cheap_finalise());
    }
    else {
#if defined(USETBB)
      {
	bool swap_and_transpose=NoSort ? 0 : (lhs.nz()+rhs.nz()<col_ptrs[rhs.cols()]);
	if (swap_and_transpose){
	  cs_free(col_ptrs); //need the swapped transpose version...
	  col_ptrs=nullptr;
	  SparseMatrix lhsTrans(lhs.copy_transpose()); //temp transpose
	  SparseMatrix rhsTrans(rhs.copy_transpose()); //temp transpose
	  col_ptrs=allocate_nonzero_totals(rhsTrans.m_array,lhsTrans.m_array);
	  reduce_sparse body(rhsTrans.m_array,lhsTrans.m_array,col_ptrs,1); //don't sort rows
	  tbb::parallel_reduce(tbb::blocked_range<SuiteSparse_long>(0,lhsTrans.m_array->n),body);
	  return std::move(SparseMatrix(body.ans_,1).transpose()); //return the transpose
	}
	else {
	  reduce_sparse body(lhs.m_array,rhs.m_array,col_ptrs,NoSort);
	  tbb::parallel_reduce(tbb::blocked_range<SuiteSparse_long>(0,rhs.m_array->n),body);
	  return SparseMatrix(body.ans_,1);
	}
      }
#else
      {
	//no TBB
	//std::cout << "NOT THREADED" <<std::endl;

	Sparseint ans_nz=col_ptrs[rhs.cols()];
	cs_free(col_ptrs);

	bool swap_and_transpose=NoSort ? 0 : (lhs.nz()+rhs.nz()<ans_nz);
	Sparseint fdim= swap_and_transpose ? rhs.cols() : lhs.rows();
	bool* Flag=new bool[fdim];
	std::fill(Flag, Flag+fdim, 0);
	const SparseMatrix* Aptr=nullptr;
	const SparseMatrix* Bptr=nullptr;
	SparseMatrix* temp1=nullptr;
	SparseMatrix* temp2=nullptr;
	std::complex<double>* W=new std::complex<double>[fdim];
	if (swap_and_transpose){ //C=(A * B)' A=rhs', B=lhs'
	  temp1=new SparseMatrix(rhs.copy_transpose());
	  temp2=new SparseMatrix(lhs.copy_transpose());
	  Aptr=temp1;
	  Bptr=temp2;
	  //delete[] Flag;
	  //Flag=new bool[rhs.cols()];///!!!!Important to reallocate Flag appropriately.
	}
	else{
	  Aptr=&lhs;
	  Bptr=&rhs;
	}
	SparseMatrix ans(Aptr->rows(),Bptr->cols(),ans_nz,1);
	ans_nz=0;
	std::fill(Flag, Flag+Aptr->rows(), 0);
	for (Sparseint k=0;k<Bptr->cols();++k){
	  ans.put_p(k)=ans_nz;
	  for (Sparseint rp=Bptr->get_p(k);rp<Bptr->get_p(k+1);++rp){
	    Sparseint j=Bptr->get_i(rp);
	    std::complex<double> Bx=Bptr->get_x(rp);
	    for (Sparseint lp=Aptr->get_p(j);lp<Aptr->get_p(j+1);++lp){
	      Sparseint i=Aptr->get_i(lp);
	      if (!Flag[i]){//don't have it yet...
		Flag[i]=1;
		ans.put_i(ans_nz++)=i;
		W[i]=Aptr->get_x(lp)*Bx;
	      }
	      else {
		W[i]+=Aptr->get_x(lp)*Bx;//scatter
	      }
	    }
	  }
	  for (Sparseint cp = ans.get_p(k) ; cp < ans_nz ; cp++){//gather
	    Sparseint i = ans.get_i(cp) ;
	    ans.put_x(cp) = W[i] ;
	    Flag[i]=0;
	  }
	}
	ans.put_p(Bptr->cols())=ans_nz;

	delete temp1;
	delete temp2;
	delete[] W;
	delete[] Flag;

	//if we want sorted by row then always need at least one transpose
	if (!NoSort) {
	  if (!swap_and_transpose){order_rows(ans.m_array);} //if we didn't do the swap transpose, then we need another now for row ordering
	  else {ans.transpose();}
	}
	return ans;
      }
#endif
    }
  }

  SparseMatrix dense_dense_multiply(const SparseMatrix& lhs, const SparseMatrix& rhs){
    SparseMatrix ans(lhs.rows(),rhs.cols(),lhs.rows()*rhs.cols(),1);
    Sparseint cptr=0;
    for (Sparseint c=0;c<rhs.cols();++c){//fill cols
      ans.put_p(c)=cptr;
      for (Sparseint i=0;i<lhs.rows();++i){
	ans.put_i(i+cptr)=i;
      }
      cptr+=lhs.rows();
    }
    ans.put_p(rhs.cols())=cptr;
    //use blas
    densefuncs::gemm(lhs.rows(), rhs.rows(), rhs.cols(), 1.0, 0.0, lhs.m_array->x, rhs.m_array->x, ans.m_array->x);
    return ans;
  }

}
