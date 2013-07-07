#ifndef VIENNACL_GENERATOR_GENERATE_HPP
#define VIENNACL_GENERATOR_GENERATE_HPP

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


/** @file viennacl/generator/generate.hpp
    @brief User interface
*/

#include <vector>
#include "viennacl/scheduler/forwards.h"
#include "viennacl/generator/detail.hpp"
namespace viennacl{

  namespace generator{

    using namespace viennacl::scheduler;

    void make_program_name(statement::container_type const & array, std::size_t index, std::string & str){
      if(array[index].op_family_==OPERATION_UNARY_TYPE_FAMILY){
        str += detail::generate(array[index].op_type_);
        str += '(';
        if(array[index].lhs_type_==COMPOSITE_OPERATION_TYPE)
          make_program_name(array, array[index].lhs_.node_index_, str);
        else
          str += detail::generate(array[index].lhs_type_);
        str += ')';
      }
      if(array[index].op_family_==OPERATION_BINARY_TYPE_FAMILY){
        str += '(';
        if(array[index].lhs_type_==COMPOSITE_OPERATION_TYPE)
          make_program_name(array, array[index].lhs_.node_index_, str);
        else
          str += detail::generate(array[index].lhs_type_);
        str += detail::generate(array[index].op_type_);
        if(array[index].rhs_type_==COMPOSITE_OPERATION_TYPE)
          make_program_name(array, array[index].rhs_.node_index_, str);
        else
          str += detail::generate(array[index].rhs_type_);
        str += ')';
      }

    }

    std::string make_program_name(scheduler::statement const & statement){
      std::string res;
      make_program_name(statement.array(), 0, res);
      return res;
    }

    template<class ExpressionsContainer>
    std::string make_program_string(ExpressionsContainer const & statements){
      return "";
    }

  }

}
#endif
