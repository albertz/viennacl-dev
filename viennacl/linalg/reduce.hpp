#ifndef VIENNACL_LINALG_REDUCE_HPP_
#define VIENNACL_LINALG_REDUCE_HPP_

/* =========================================================================
   Copyright (c) 2010-2013, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/reduce.hpp
    @brief Generic interface for the computation of inner products. See viennacl/linalg/vector_operations.hpp for implementations.
*/

#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/meta/tag_of.hpp"
#include "viennacl/meta/result_of.hpp"

namespace viennacl
{

  namespace linalg
  {

    // ----------------------------------------------------
    template<typename OP, typename NumericT>
    viennacl::scalar_expression< const vector_base<NumericT>, const vector_base<NumericT>, viennacl::op_reduce<OP> >
    reduce(vector_base<NumericT> const & vector)
    {
      //std::cout << "viennacl .. " << std::endl;
      return viennacl::scalar_expression< const vector_base<NumericT>,
                                          const vector_base<NumericT>,
                                          viennacl::op_reduce<OP> >(vector, vector);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif


