#ifndef VIENNACL_GENERATOR_DETAIL_HPP
#define VIENNACL_GENERATOR_DETAIL_HPP

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


/** @file viennacl/generator/detail.hpp
    @brief Internal helper for generating strings
*/

#include "viennacl/scheduler/forwards.h"

namespace viennacl{

  namespace generator{

    namespace detail{

        using namespace viennacl::scheduler;


#define MAKE_CASE(ref,str) if(arg==ref) return str;

        const char * generate(operation_node_type arg){
            // unary expression
            MAKE_CASE(OPERATION_UNARY_ABS_TYPE, "abs")

            // binary expression
            MAKE_CASE(OPERATION_BINARY_ASSIGN_TYPE, "=")
            MAKE_CASE(OPERATION_BINARY_ADD_TYPE, "+")

            throw "missing operator";
        }

        const char * generate(statement_node_type arg){
            MAKE_CASE(COMPOSITE_OPERATION_TYPE,"")
            //vector
            MAKE_CASE(VECTOR_FLOAT_TYPE,"vecf")

            throw "missing node";
        }

#undef MAKE_CASE

        template<class Fun>
        void for_each(viennacl::scheduler::statement const & statements, Fun const & fun){
            for(typename viennacl::scheduler::statement::container_type::const_iterator it = statements.array().begin() ; it != statements.array().end() ; ++it){
              fun(it->lhs_type_);
              fun(it->op_type_);
            }
            fun(statements.array().back().rhs_type_);
        }

    }

  }

}
#endif
