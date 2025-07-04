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


//! \addtogroup op_misc
//! @{



template<typename T1>
inline
void
op_real::apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_real>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X.m);
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
    
  out.set_size(n_rows, n_cols);
  
  T* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = std::real( A[i] );
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      *out_mem = std::real( P.at(row,col) );
      out_mem++;
      }
    }
  }



template<typename T1>
inline
void
op_real::apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_real>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const ProxyCube<T1> P(X.m);
  
  const uword n_rows   = P.get_n_rows();
  const uword n_cols   = P.get_n_cols();
  const uword n_slices = P.get_n_slices();
    
  out.set_size(n_rows, n_cols, n_slices);
  
  T* out_mem = out.memptr();

  if(ProxyCube<T1>::use_at == false)
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = std::real( A[i] );
      }
    }
  else
    {
    for(uword slice=0; slice < n_slices; ++slice)
    for(uword col=0;   col   < n_cols;   ++col  )
    for(uword row=0;   row   < n_rows;   ++row  )
      {
      *out_mem = std::real( P.at(row,col,slice) );
      out_mem++;
      }
    }
  }



template<typename T1>
inline
void
op_imag::apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_imag>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> P(X.m);
  
  if(is_cx<eT>::no)  { out.zeros(P.get_n_rows(), P.get_n_cols()); return; }
  
  // aliasing not possible at this point, as eT must be std::complex
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
    
  out.set_size(n_rows, n_cols);
  
  T* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = access::tmp_imag( A[i] );
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      *out_mem = access::tmp_imag( P.at(row,col) );
      out_mem++;
      }
    }
  }



template<typename T1>
inline
void
op_imag::apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_imag>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const ProxyCube<T1> P(X.m);
  
  if(is_cx<eT>::no)  { out.zeros(P.get_n_rows(), P.get_n_cols(), P.get_n_slices()); return; }
  
  // aliasing not possible at this point, as eT must be std::complex
  
  const uword n_rows   = P.get_n_rows();
  const uword n_cols   = P.get_n_cols();
  const uword n_slices = P.get_n_slices();
    
  out.set_size(n_rows, n_cols, n_slices);
  
  T* out_mem = out.memptr();

  if(ProxyCube<T1>::use_at == false)
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = access::tmp_imag( A[i] );
      }
    }
  else
    {
    for(uword slice=0; slice < n_slices; ++slice)
    for(uword col=0;   col   < n_cols;   ++col  )
    for(uword row=0;   row   < n_rows;   ++row  )
      {
      *out_mem = access::tmp_imag( P.at(row,col,slice) );
      out_mem++;
      }
    }
  }



template<typename T1>
inline
void
op_abs::apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_abs>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X.m);
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
  T* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i < n_elem; ++i)
        {
        out_mem[i] = std::abs( A[i] );
        }
      }
    #else
      {
      for(uword i=0; i < n_elem; ++i)
        {
        out_mem[i] = std::abs( A[i] );
        }
      }
    #endif
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      *out_mem = std::abs( P.at(row,col) );
      out_mem++;
      }
    }
  }



template<typename T1>
inline
void
op_abs::apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_abs>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const ProxyCube<T1> P(X.m);
  
  const uword n_rows   = P.get_n_rows();
  const uword n_cols   = P.get_n_cols();
  const uword n_slices = P.get_n_slices();
    
  out.set_size(n_rows, n_cols, n_slices);
  
  T* out_mem = out.memptr();

  if(ProxyCube<T1>::use_at == false)
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i < n_elem; ++i)
        {
        out_mem[i] = std::abs( A[i] );
        }
      }
    #else
      {
      for(uword i=0; i < n_elem; ++i)
        {
        out_mem[i] = std::abs( A[i] );
        }
      }
    #endif
    }
  else
    {
    for(uword slice=0; slice < n_slices; ++slice)
    for(uword col=0;   col   < n_cols;   ++col  )
    for(uword row=0;   row   < n_rows;   ++row  )
      {
      *out_mem = std::abs( P.at(row,col,slice) );
      out_mem++;
      }
    }
  }



template<typename T1>
inline
void
op_arg::apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_arg>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> P(X.m);
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
  T* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = arma_arg<eT>::eval( A[i] );
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      *out_mem = arma_arg<eT>::eval( P.at(row,col) );
      out_mem++;
      }
    }
  }



template<typename T1>
inline
void
op_arg::apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_arg>& X )
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const ProxyCube<T1> P(X.m);
  
  const uword n_rows   = P.get_n_rows();
  const uword n_cols   = P.get_n_cols();
  const uword n_slices = P.get_n_slices();
  
  out.set_size(n_rows, n_cols, n_slices);
  
  T* out_mem = out.memptr();

  if(ProxyCube<T1>::use_at == false)
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const uword   n_elem  = P.get_n_elem();
          ea_type A       = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = arma_arg<eT>::eval( A[i] );
      }
    }
  else
    {
    for(uword slice=0; slice < n_slices; ++slice)
    for(uword col=0;   col   < n_cols;   ++col  )
    for(uword row=0;   row   < n_rows;   ++row  )
      {
      *out_mem = arma_arg<eT>::eval( P.at(row,col,slice) );
      out_mem++;
      }
    }
  }



template<typename eT, typename T1>
inline
void
op_replace::apply(Mat<eT>& out, const mtOp<eT,T1,op_replace>& in)
  {
  arma_debug_sigprint();
  
  const eT old_val = in.aux;
  const eT new_val = in.aux_out_eT;
  
  out = in.m;
  
  out.replace(old_val, new_val);
  }



template<typename eT, typename T1>
inline
void
op_replace::apply(Cube<eT>& out, const mtOpCube<eT,T1,op_replace>& in)
  {
  arma_debug_sigprint();
  
  const eT old_val = in.aux;
  const eT new_val = in.aux_out_eT;
  
  out = in.m;
  
  out.replace(old_val, new_val);
  }



//! @}
