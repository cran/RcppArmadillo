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


//! \addtogroup mtOp
//! @{


struct mtOp_dual_aux_indicator {};


template<typename out_eT, typename T1, typename op_type>
class mtOp : public Base< out_eT, mtOp<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  typedef typename T1::elem_type                in_eT;
  
  static constexpr bool is_row  = op_type::template traits<T1>::is_row;
  static constexpr bool is_col  = op_type::template traits<T1>::is_col;
  static constexpr bool is_xvec = op_type::template traits<T1>::is_xvec;
  
  inline explicit mtOp(const T1& in_m);
  inline          mtOp(const T1& in_m, const in_eT in_aux);
  inline          mtOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
  inline          mtOp(const T1& in_m, const in_eT in_aux,         const uword in_aux_uword_a, const uword in_aux_uword_b);
  
  inline          mtOp(const char junk, const T1& in_m, const out_eT in_aux);
  
  inline          mtOp(const mtOp_dual_aux_indicator&, const T1& in_m, const in_eT in_aux_a, const out_eT in_aux_b);
  
  inline         ~mtOp();
  
  template<typename eT2>
  inline bool is_alias(const Mat<eT2>& X) const;
  
  const T1&    m;            //!< the operand; must be derived from Base
        in_eT  aux;          //!< auxiliary data, using the element type as used by T1
        out_eT aux_out_eT;   //!< auxiliary data, using the element type as specified by the out_eT template parameter
        uword  aux_uword_a;  //!< auxiliary data, uword format
        uword  aux_uword_b;  //!< auxiliary data, uword format
  };



//! @}
