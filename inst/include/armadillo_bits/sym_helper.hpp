// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup sym_helper
//! @{


namespace sym_helper
{

// computationally inexpensive algorithm to guess whether a matrix is positive definite:
// (1) ensure the matrix is symmetric/hermitian (within a tolerance)
// (2) ensure the diagonal entries are real and greater than zero
// (3) ensure that the value with largest modulus is on the main diagonal
// (4) ensure rudimentary diagonal dominance: (real(A_ii) + real(A_jj)) > 2*abs(real(A_ij))
// the above conditions are necessary, but not sufficient;
// doing it properly would be too computationally expensive for our purposes
// more info:
// http://mathworld.wolfram.com/PositiveDefiniteMatrix.html
// http://mathworld.wolfram.com/DiagonallyDominantMatrix.html
  
template<typename eT>
inline
typename enable_if2<is_cx<eT>::no, bool>::result
guess_sympd_worker(const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  // NOTE: assuming A is square-sized
  
  const eT tol = eT(100) * std::numeric_limits<eT>::epsilon();  // allow some leeway
  
  const uword N = A.n_rows;
  
  const eT* A_mem = A.memptr();
  const eT* A_col = A_mem;
  
  bool diag_below_tol = true;
  
  eT max_diag = eT(0);
  
  for(uword j=0; j < N; ++j)
    {
    const eT A_jj = A_col[j];
    
    if(        A_jj  <= eT(0))  { return false; }
    if(arma_isnonfinite(A_jj))  { return false; }
    
    if(A_jj >= tol)  { diag_below_tol = false; }
    
    max_diag = (A_jj > max_diag) ? A_jj : max_diag;
    
    A_col += N;
    }
  
  if(diag_below_tol)  { return false; }  // assume matrix is suspect if all diagonal elements are close to zero
  
  A_col = A_mem;
  
  const uword Nm1 = N-1;
  const uword Np1 = N+1;
  
  for(uword j=0; j < Nm1; ++j)
    {
    const eT A_jj = A_col[j];
    
    const uword jp1      = j+1;
    const eT*   A_ji_ptr = &(A_mem[j   + jp1*N]);  // &(A.at(j,jp1));
    const eT*   A_ii_ptr = &(A_mem[jp1 + jp1*N]);
    
    for(uword i=jp1; i < N; ++i)
      {
      const eT A_ij = A_col[i];
      const eT A_ji = (*A_ji_ptr);
      
      const eT A_ij_abs = (std::abs)(A_ij);
      const eT A_ji_abs = (std::abs)(A_ji);
      
      // if( (A_ij_abs >= max_diag) || (A_ji_abs >= max_diag) )  { return false; }
      if(A_ij_abs >= max_diag)  { return false; }
      
      const eT A_delta   = (std::abs)(A_ij - A_ji);
      const eT A_abs_max = (std::max)(A_ij_abs, A_ji_abs);
      
      if( (A_delta > tol) && (A_delta > (A_abs_max*tol)) )  { return false; }
      
      const eT A_ii = (*A_ii_ptr);
      
      if( (A_ij_abs + A_ij_abs) >= (A_ii + A_jj) )  { return false; }
      
      A_ji_ptr += N;
      A_ii_ptr += Np1;
      }
    
    A_col += N;
    }
  
  return true;
  }



template<typename eT>
inline
typename enable_if2<is_cx<eT>::yes, bool>::result
guess_sympd_worker(const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  // NOTE: assuming A is square-sized
  
  // NOTE: the function name is required for overloading, but is a misnomer: it processes complex hermitian matrices
  
  typedef typename get_pod_type<eT>::result T;
  
  const T tol = T(100) * std::numeric_limits<T>::epsilon();  // allow some leeway
  
  const uword N = A.n_rows;
  
  const eT* A_mem = A.memptr();
  const eT* A_col = A_mem;
  
  bool diag_below_tol = true;
  
  T max_diag = T(0);
  
  for(uword j=0; j < N; ++j)
    {
    const eT& A_jj      = A_col[j];
    const  T  A_jj_r    = std::real(A_jj  );
    const  T  A_jj_i    = std::imag(A_jj  );
    const  T  A_jj_rabs =  std::abs(A_jj_r);
    const  T  A_jj_iabs =  std::abs(A_jj_i);
    
    if(        A_jj_r  <= T(0) )  { return false; }  // real should be positive
    if(arma_isnonfinite(A_jj_r))  { return false; }
    
    if(A_jj_iabs > tol      )  { return false; }  // imag should be approx zero
    if(A_jj_iabs > A_jj_rabs)  { return false; }  // corner case: real and imag are close to zero, and imag is dominant
    
    if(A_jj_r >= tol)  { diag_below_tol = false; }
    
    max_diag = (A_jj_r > max_diag) ? A_jj_r : max_diag;
    
    A_col += N;
    }
  
  if(diag_below_tol)  { return false; }  // assume matrix is suspect if all diagonal elements are close to zero
  
  const T square_max_diag = max_diag * max_diag;
  
  if(arma_isnonfinite(square_max_diag))  { return false; }
  
  A_col = A_mem;
  
  const uword Nm1 = N-1;
  const uword Np1 = N+1;
  
  for(uword j=0; j < Nm1; ++j)
    {
    const uword jp1       = j+1;
    const eT*   A_ji_ptr = &(A_mem[j   + jp1*N]);  // &(A.at(j,jp1));
    const eT*   A_ii_ptr = &(A_mem[jp1 + jp1*N]);
    
    const T A_jj_real = std::real(A_col[j]);
    
    for(uword i=jp1; i < N; ++i)
      {
      const eT& A_ij      = A_col[i];
      const  T  A_ij_real = std::real(A_ij);
      const  T  A_ij_imag = std::imag(A_ij);
      
      // avoid using std::abs(), as that is time consuming due to division and std::sqrt()
      const T square_A_ij_abs = (A_ij_real * A_ij_real) + (A_ij_imag * A_ij_imag);
      
      if(arma_isnonfinite(square_A_ij_abs))  { return false; }
      
      if(square_A_ij_abs >= square_max_diag)  { return false; }
      
      const T A_ij_real_abs = (std::abs)(A_ij_real);
      const T A_ij_imag_abs = (std::abs)(A_ij_imag);
      
      
      const eT& A_ji      = (*A_ji_ptr);
      const  T  A_ji_real = std::real(A_ji);
      const  T  A_ji_imag = std::imag(A_ji);
      
      const T A_ji_real_abs = (std::abs)(A_ji_real);
      const T A_ji_imag_abs = (std::abs)(A_ji_imag);
      
      const T A_real_delta   = (std::abs)(A_ij_real - A_ji_real);
      const T A_real_abs_max = (std::max)(A_ij_real_abs, A_ji_real_abs);
      
      if( (A_real_delta > tol) && (A_real_delta > (A_real_abs_max*tol)) )  { return false; }
      
      
      const T A_imag_delta   = (std::abs)(A_ij_imag + A_ji_imag);  // take into account complex conjugate
      const T A_imag_abs_max = (std::max)(A_ij_imag_abs, A_ji_imag_abs);
      
      if( (A_imag_delta > tol) && (A_imag_delta > (A_imag_abs_max*tol)) )  { return false; }
      
      
      const T A_ii_real = std::real(*A_ii_ptr);
      
      if( (A_ij_real_abs + A_ij_real_abs) >= (A_ii_real + A_jj_real) )  { return false; }
      
      A_ji_ptr += N;
      A_ii_ptr += Np1;
      }
    
    A_col += N;
    }
  
  return true;
  }



template<typename eT>
inline
bool
guess_sympd(const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  // analyse matrices with size >= 4x4
  
  if((A.n_rows != A.n_cols) || (A.n_rows < uword(4)))  { return false; }
  
  return guess_sympd_worker(A);
  }



template<typename eT>
inline
bool
guess_sympd(const Mat<eT>& A, const uword min_n_rows)
  {
  arma_debug_sigprint();
  
  if((A.n_rows != A.n_cols) || (A.n_rows < min_n_rows))  { return false; }
  
  return guess_sympd_worker(A);
  }



//



template<typename eT>
inline
typename enable_if2<is_cx<eT>::no, bool>::result
is_approx_sym_worker(const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  const eT tol = eT(100) * std::numeric_limits<eT>::epsilon();  // allow some leeway
  
  const uword N = A.n_rows;
  
  const eT* A_mem = A.memptr();
  const eT* A_col = A_mem;
  
  bool diag_below_tol = true;
  
  for(uword j=0; j < N; ++j)
    {
    const eT A_jj = A_col[j];
    
    if(arma_isnonfinite(A_jj))  { return false; }
    
    if(std::abs(A_jj) >= tol)  { diag_below_tol = false; }
    
    A_col += N;
    }
  
  if(diag_below_tol)  { return false; }  // assume matrix is suspect if all diagonal elements are close to zero
  
  A_col = A_mem;
  
  const uword Nm1 = N-1;
  
  for(uword j=0; j < Nm1; ++j)
    {
    const uword jp1      = j+1;
    const eT*   A_ji_ptr = &(A_mem[j + jp1*N]);  // &(A.at(j,jp1));
    
    for(uword i=jp1; i < N; ++i)
      {
      const eT A_ij = A_col[i];
      const eT A_ji = (*A_ji_ptr);
      
      const eT A_ij_abs = (std::abs)(A_ij);
      const eT A_ji_abs = (std::abs)(A_ji);
      
      const eT A_delta   = (std::abs)(A_ij - A_ji);
      const eT A_abs_max = (std::max)(A_ij_abs, A_ji_abs);
      
      if( (A_delta > tol) && (A_delta > (A_abs_max*tol)) )  { return false; }
      
      A_ji_ptr += N;
      }
    
    A_col += N;
    }
  
  return true;
  }



template<typename eT>
inline
typename enable_if2<is_cx<eT>::yes, bool>::result
is_approx_sym_worker(const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  // NOTE: the function name is required for overloading, but is a misnomer: it processes complex hermitian matrices
  
  typedef typename get_pod_type<eT>::result T;
  
  const T tol = T(100) * std::numeric_limits<T>::epsilon();  // allow some leeway
  
  const uword N = A.n_rows;
  
  const eT* A_mem = A.memptr();
  const eT* A_col = A_mem;
  
  bool diag_below_tol = true;
  
  // ensure diagonal has approx real-only elements
  for(uword j=0; j < N; ++j)
    {
    const eT& A_jj      = A_col[j];
    const  T  A_jj_r    = std::real(A_jj  );
    const  T  A_jj_i    = std::imag(A_jj  );
    const  T  A_jj_rabs =  std::abs(A_jj_r);
    const  T  A_jj_iabs =  std::abs(A_jj_i);
    
    if(A_jj_iabs > tol      )  { return false; }  // imag should be approx zero
    if(A_jj_iabs > A_jj_rabs)  { return false; }  // corner case: real and imag are close to zero, and imag is dominant
    
    if(arma_isnonfinite(A_jj_r))  { return false; }
    
    if(A_jj_rabs >= tol)  { diag_below_tol = false; }
    
    A_col += N;
    }
  
  if(diag_below_tol)  { return false; }  // assume matrix is suspect if all diagonal elements are close to zero
  
  A_col = A_mem;
  
  const uword Nm1 = N-1;
  
  for(uword j=0; j < Nm1; ++j)
    {
    const uword jp1      = j+1;
    const eT*   A_ji_ptr = &(A_mem[j + jp1*N]);  // &(A.at(j,jp1));
    
    for(uword i=jp1; i < N; ++i)
      {
      const eT& A_ij      = A_col[i];
      const  T  A_ij_real = std::real(A_ij);
      const  T  A_ij_imag = std::imag(A_ij);
      
      const T A_ij_real_abs = (std::abs)(A_ij_real);
      const T A_ij_imag_abs = (std::abs)(A_ij_imag);
      
      const eT& A_ji      = (*A_ji_ptr);
      const  T  A_ji_real = std::real(A_ji);
      const  T  A_ji_imag = std::imag(A_ji);
      
      const T A_ji_real_abs = (std::abs)(A_ji_real);
      const T A_ji_imag_abs = (std::abs)(A_ji_imag);
      
      const T A_real_delta   = (std::abs)(A_ij_real - A_ji_real);
      const T A_real_abs_max = (std::max)(A_ij_real_abs, A_ji_real_abs);
      
      if( (A_real_delta > tol) && (A_real_delta > (A_real_abs_max*tol)) )  { return false; }
      
      const T A_imag_delta   = (std::abs)(A_ij_imag + A_ji_imag);  // take into account complex conjugate
      const T A_imag_abs_max = (std::max)(A_ij_imag_abs, A_ji_imag_abs);
      
      if( (A_imag_delta > tol) && (A_imag_delta > (A_imag_abs_max*tol)) )  { return false; }
      
      A_ji_ptr += N;
      }
    
    A_col += N;
    }
  
  return true;
  }



template<typename eT>
inline
bool
is_approx_sym(const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  // analyse matrices with size >= 4x4
  
  if((A.n_rows != A.n_cols) || (A.n_rows < uword(4)))  { return false; }
  
  return is_approx_sym_worker(A);
  }



template<typename eT>
inline
bool
is_approx_sym(const Mat<eT>& A, const uword min_n_rows)
  {
  arma_debug_sigprint();
  
  if((A.n_rows != A.n_cols) || (A.n_rows < min_n_rows))  { return false; }
  
  return is_approx_sym_worker(A);
  }



//



template<typename eT>
inline
bool
check_diag_imag(const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  // NOTE: assuming matrix A is square-sized
  
  typedef typename get_pod_type<eT>::result T;
  
  const T tol = T(10000) * std::numeric_limits<T>::epsilon();  // allow some leeway
  
  const eT* colmem = A.memptr();
  
  const uword N = A.n_rows;
  
  for(uword i=0; i<N; ++i)
    {
    const eT& A_ii      = colmem[i];
    const  T  A_ii_imag = access::tmp_imag(A_ii);
    
    if(std::abs(A_ii_imag) > tol)  { return false; }
    
    colmem += N;
    }
  
  return true;
  }



}  // end of namespace sym_helper


//! @}
