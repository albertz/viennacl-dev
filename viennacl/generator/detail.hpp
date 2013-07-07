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
    @brief Internal upper layer to the scheduler
*/

#include "CL/cl.h"

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/utils.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{

      using namespace viennacl::scheduler;

      cl_mem get_handle(statement_node_type type, lhs_rhs_element element){
#define MAKE_CASE(ref,unmap_type, raw_pointer) if(type==ref) return static_cast<unmap_type *>(element.raw_pointer)->handle().opencl_handle();
        //vector
        MAKE_CASE(VECTOR_FLOAT_TYPE, vector_base<float>, vector_float_);

        throw "unrecognized type";
#undef MAKE_CASE
      }

      std::string generate_scalartype(statement_node_type type){
#define MAKE_CASE(ref, scalartype) if(type==ref) return scalartype;
        //vector
        MAKE_CASE(VECTOR_FLOAT_TYPE, "float");

        throw "unrecognized type";
#undef MAKE_CASE
      }

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

      template<class TraversalFunctor>
      void traverse(statement::container_type const & array, std::size_t index, TraversalFunctor const & fun){
        statement::value_type const & element = array[index];
        if(element.op_family_==OPERATION_UNARY_TYPE_FAMILY){
          fun.call_on_op(element.op_family_, element.op_type_);
          fun.call_before_expansion();
          if(element.lhs_type_==COMPOSITE_OPERATION_TYPE)
            traverse(array, element.lhs_.node_index_, fun);
          else
            fun.call_on_leaf(element.lhs_type_family_, element.lhs_type_, element.lhs_);
          fun.call_after_expansion();
        }
        if(element.op_family_==OPERATION_BINARY_TYPE_FAMILY){
          fun.call_before_expansion();
          if(element.lhs_type_==COMPOSITE_OPERATION_TYPE)
            traverse(array, element.lhs_.node_index_, fun);
          else
            fun.call_on_leaf(element.lhs_type_family_, element.lhs_type_, element.lhs_);
          fun.call_on_op(element.op_family_, element.op_type_);
          if(element.rhs_type_==COMPOSITE_OPERATION_TYPE)
            traverse(array, element.rhs_.node_index_, fun);
          else
            fun.call_on_leaf(element.rhs_type_family_, element.rhs_type_, element.rhs_);
          fun.call_after_expansion();
        }
      }

      class generation_traversal{
          std::string & str_;
        public:
          generation_traversal(std::string & str) : str_(str) { }
          void call_on_leaf(statement_node_type_family, statement_node_type type, lhs_rhs_element) const { str_ += detail::generate(type); }
          void call_on_op(operation_node_type_family, operation_node_type type) const { str_ += detail::generate(type); }
          void call_before_expansion() const { str_ += '('; }
          void call_after_expansion() const { str_ += ')'; }
      };

      class header_generation_traversal{
          mutable unsigned int current_;
          std::map<cl_mem, unsigned int> & memory_;
          std::string & str_;
          void header_value_generation(statement_node_type_family, statement_node_type type, lhs_rhs_element) const{
            str_ += detail::generate_scalartype(type) + ' ' + "arg" + utils::to_string(current_++) + ",";
          }
          void header_pointer_generation(statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element) const {
            if(memory_.insert(std::make_pair(detail::get_handle(type, element), current_)).second)
              str_ += "__global " + detail::generate_scalartype(type) + "* arg" + utils::to_string(current_++) + ",";
          }
        public:
          header_generation_traversal(std::map<cl_mem, unsigned int> & memory, std::string & str) : current_(0), memory_(memory), str_(str) { }
          void call_on_leaf(statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element) const {
            if(type_family==HOST_SCALAR_TYPE_FAMILY)
              header_value_generation(type_family,type,element);
            else
              header_pointer_generation(type_family,type,element);
          }
          void call_on_op(operation_node_type_family, operation_node_type type) const {  }
          void call_before_expansion() const { }
          void call_after_expansion() const { }
      };


    }

  }

}
#endif
