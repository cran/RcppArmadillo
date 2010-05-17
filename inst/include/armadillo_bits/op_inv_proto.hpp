// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_inv
//! @{



//! 'invert matrix' operation

class op_inv
  {
  public:
  
  // mat

  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv>& in);

  template<typename T1>
  inline static void apply_diag(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& X);

  };



//! @}
