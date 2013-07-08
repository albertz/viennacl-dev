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
            fun.call_on_leaf(index, element.lhs_type_family_, element.lhs_type_, element.lhs_);
          fun.call_after_expansion();
        }
        if(element.op_family_==OPERATION_BINARY_TYPE_FAMILY){
          fun.call_before_expansion();
          if(element.lhs_type_==COMPOSITE_OPERATION_TYPE)
            traverse(array, element.lhs_.node_index_, fun);
          else
            fun.call_on_leaf(index, element.lhs_type_family_, element.lhs_type_, element.lhs_);
          fun.call_on_op(element.op_family_, element.op_type_);
          if(element.rhs_type_==COMPOSITE_OPERATION_TYPE)
            traverse(array, element.rhs_.node_index_, fun);
          else
            fun.call_on_leaf(index, element.rhs_type_family_, element.rhs_type_, element.rhs_);
          fun.call_after_expansion();
        }
      }

      class traversal_functor{
        public:
          void call_on_op(operation_node_type_family, operation_node_type type) const { }
          void call_before_expansion() const { }
          void call_after_expansion() const { }
      };

      struct symbolic_container{
          std::string scalartype_;
          std::string name_;
      };

      typedef std::map<std::size_t, detail::symbolic_container> mapping_type;

      template<class Fun>
      class expression_generation_traversal : public traversal_functor{
          utils::kernel_generation_stream & stream_;
          mapping_type const & mapping_;
          Fun fun_;
        public:
          expression_generation_traversal(utils::kernel_generation_stream & stream, mapping_type const & mapping, Fun const & fun) : stream_(stream), mapping_(mapping), fun_(fun){ }
          void call_on_leaf(std::size_t index, statement_node_type_family, statement_node_type type, lhs_rhs_element) const { stream_ << fun_(mapping_.at(index)); }
          void call_on_op(operation_node_type_family, operation_node_type type) const { stream_ << detail::generate(type); }
          void call_before_expansion() const { stream_ << '('; }
          void call_after_expansion() const { stream_ << ')'; }
      };

      class add_suffix_to_name{
          std::string suffix_;
        public:
          add_suffix_to_name(std::string const & suffix) : suffix_(suffix){ }
          std::string operator()(detail::symbolic_container const & sym) const { return sym.name_ + suffix_; }
      };

    }

  }

}
#endif
