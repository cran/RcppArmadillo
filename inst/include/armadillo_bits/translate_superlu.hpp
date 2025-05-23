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



#if defined(ARMA_USE_SUPERLU)

//! \namespace superlu namespace for SuperLU functions
namespace superlu
  {
  
  template<typename eT>
  inline
  void
  gssv(superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, superlu::int_t* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    if(is_float<eT>::value)
      {
      arma_wrapper(sgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    else
    if(is_double<eT>::value)
      {
      arma_wrapper(dgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    else
    if(is_cx_float<eT>::value)
      {
      arma_wrapper(cgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    else
    if(is_cx_double<eT>::value)
      {
      arma_wrapper(zgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  gssvx(
        superlu_options_t* opts,
        SuperMatrix* A,
        int* perm_c, int* perm_r,
        int* etree, char* equed,
        typename get_pod_type<eT>::result* R, typename get_pod_type<eT>::result* C,
        SuperMatrix* L, SuperMatrix* U,
        void* work, superlu::int_t lwork,
        SuperMatrix* B, SuperMatrix* X,
        typename get_pod_type<eT>::result* rpg, typename get_pod_type<eT>::result* rcond,
        typename get_pod_type<eT>::result* ferr, typename get_pod_type<eT>::result* berr,
        GlobalLU_t* glu, mem_usage_t* mu, SuperLUStat_t* stat, superlu::int_t* info
       )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(sgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(dgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    else
    if(is_cx_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(cgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    else
    if(is_cx_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(zgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  gstrf(superlu_options_t* options,
        SuperMatrix* A,
        int relax,
        int panel_size, int* etree,
        void* work, superlu::int_t lwork,
        int* perm_c, int* perm_r,
        SuperMatrix* L, SuperMatrix* U,
        GlobalLU_t* Glu, SuperLUStat_t* stat, superlu::int_t* info
       )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));

    if(is_float<eT>::value)
      {
      arma_wrapper(sgstrf)(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, Glu, stat, info);
      }
    else
    if(is_double<eT>::value)
      {
      arma_wrapper(dgstrf)(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, Glu, stat, info);
      }
    else
    if(is_cx_float<eT>::value)
      {
      arma_wrapper(cgstrf)(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, Glu, stat, info);
      }
    else
    if(is_cx_double<eT>::value)
      {
      arma_wrapper(zgstrf)(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, Glu, stat, info);
      }
    }



  template<typename eT>
  inline
  void
  gstrs(trans_t trans,
        SuperMatrix* L, SuperMatrix* U,
        int* perm_c, int* perm_r,
        SuperMatrix* B, SuperLUStat_t* stat, int* info
       )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));

    if(is_float<eT>::value)
      {
      arma_wrapper(sgstrs)(trans, L, U, perm_c, perm_r, B, stat, info);
      }
    else
    if(is_double<eT>::value)
      {
      arma_wrapper(dgstrs)(trans, L, U, perm_c, perm_r, B, stat, info);
      }
    else
    if(is_cx_float<eT>::value)
      {
      arma_wrapper(cgstrs)(trans, L, U, perm_c, perm_r, B, stat, info);
      }
    else
    if(is_cx_double<eT>::value)
      {
      arma_wrapper(zgstrs)(trans, L, U, perm_c, perm_r, B, stat, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  typename get_pod_type<eT>::result
  langs(char* norm, superlu::SuperMatrix* A)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    typedef typename get_pod_type<eT>::result T;
    
    if(is_float<eT>::value)
      {
      return arma_wrapper(slangs)(norm, A);
      }
    else
    if(is_double<eT>::value)
      {
      return arma_wrapper(dlangs)(norm, A);
      }
    else
    if(is_cx_float<eT>::value)
      {
      return arma_wrapper(clangs)(norm, A);
      }
    else
    if(is_cx_double<eT>::value)
      {
      return arma_wrapper(zlangs)(norm, A);
      }
    
    return T(0);  // to avoid false warnings from the compiler
    }
  
  
  
  template<typename eT>
  inline
  void
  gscon(char* norm, superlu::SuperMatrix* L, superlu::SuperMatrix* U, typename get_pod_type<eT>::result anorm, typename get_pod_type<eT>::result* rcond, superlu::SuperLUStat_t* stat, int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));

    if(is_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(sgscon)(norm, L, U, (T)anorm, (T*)rcond, stat, info);
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(dgscon)(norm, L, U, (T)anorm, (T*)rcond, stat, info);
      }
    else
    if(is_cx_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(cgscon)(norm, L, U, (T)anorm, (T*)rcond, stat, info);
      }
    else
    if(is_cx_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(zgscon)(norm, L, U, (T)anorm, (T*)rcond, stat, info);
      }
    }
  
  
  
  inline
  void
  init_stat(SuperLUStat_t* stat)
    {
    arma_wrapper(StatInit)(stat);
    }


  inline
  void
  free_stat(SuperLUStat_t* stat)
    {
    arma_wrapper(StatFree)(stat);
    }
  
  
  
  inline
  void
  set_default_opts(superlu_options_t* opts)
    {
    arma_wrapper(set_default_options)(opts);
    }
  
  
  inline
  void
  get_permutation_c(int ispec, SuperMatrix* A, int* perm_c)
    {
    arma_wrapper(get_perm_c)(ispec, A, perm_c);
    }



  inline
  void
  sp_preorder_mat(superlu_options_t* opts, SuperMatrix* A, int* perm_c, int* etree, SuperMatrix* AC)
    {
    arma_wrapper(sp_preorder)(opts, A, perm_c, etree, AC);
    }



  inline
  int
  sp_ispec_environ(int ispec)
    {
    return arma_wrapper(sp_ienv)(ispec);
    }



  inline
  void
  destroy_supernode_mat(SuperMatrix* a)
    {
    arma_wrapper(Destroy_SuperNode_Matrix)(a);
    }



  inline
  void
  destroy_compcol_mat(SuperMatrix* a)
    {
    arma_wrapper(Destroy_CompCol_Matrix)(a);
    }



  inline
  void
  destroy_compcolperm_mat(SuperMatrix* a)
    {
    arma_wrapper(Destroy_CompCol_Permuted)(a);
    }



  inline
  void
  destroy_dense_mat(SuperMatrix* a)
    {
    arma_wrapper(Destroy_SuperMatrix_Store)(a);
    }
  
  
  
  inline
  void*
  malloc(size_t N)
    {
    return arma_wrapper(superlu_malloc)(N);
    }
  
  
  
  inline
  void
  free(void* mem)
    {
    arma_wrapper(superlu_free)(mem);
    }
  
  } // namespace superlu

#endif
