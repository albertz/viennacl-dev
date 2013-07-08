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

      struct symbolic_container{
          std::string scalartype_;
          std::string name_;
      };

      class fetching_traversal{
          std::string i_;
          utils::kernel_generation_stream & stream_;
          std::vector<symbolic_container> const & symbolic_mapping_;
        public:
          fetching_traversal(utils::kernel_generation_stream & stream, std::vector<symbolic_container> const & symbolic_mapping, std::string const & i) : i_(i), stream_(stream), symbolic_mapping_(symbolic_mapping){ }

          void call_on_leaf(std::size_t index, statement_node_type_family, statement_node_type type, lhs_rhs_element) const {
            stream_ << symbolic_mapping_[index].scalartype_ << ' ' << symbolic_mapping_[index].name_ << ';' << std::endl;
          }

          void call_on_op(operation_node_type_family, operation_node_type type) const { }

          void call_before_expansion() const { }

          void call_after_expansion() const { }
      };

      class name_generation_traversal{
          std::string & str_;
        public:
          name_generation_traversal(std::string & str) : str_(str) { }
          void call_on_leaf(std::size_t, statement_node_type_family, statement_node_type type, lhs_rhs_element) const { str_ += detail::generate(type); }
          void call_on_op(operation_node_type_family, operation_node_type type) const { str_ += detail::generate(type); }
          void call_before_expansion() const { str_ += '('; }
          void call_after_expansion() const { str_ += ')'; }
      };


      class header_generation_traversal{

          mutable unsigned int current_arg_;
          std::map<cl_mem, std::size_t> & memory_;
          std::vector<symbolic_container> & symbolic_mapping_;
          std::string & str_;


          void header_value_generation(statement_node_type_family, statement_node_type type, lhs_rhs_element, symbolic_container & sym) const{
            sym.scalartype_ = detail::generate_scalartype(type);
            sym.name_ = "arg" + utils::to_string(current_arg_++);
            str_ += sym.scalartype_ + ' '  + sym.name_ + ",";
          }

          void header_pointer_generation(statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element, symbolic_container & sym) const {
            sym.scalartype_ = detail::generate_scalartype(type);
            if(memory_.insert(std::make_pair(detail::get_handle(type, element), current_arg_)).second){
              sym.name_ =  "arg" + utils::to_string(current_arg_++);
              str_ += "__global " +  sym.scalartype_ + "* " + sym.name_ + ",";
            }
            sym.name_ = "arg" + utils::to_string(memory_.at(detail::get_handle(type, element)));
          }

        public:
          header_generation_traversal(std::map<cl_mem, std::size_t> & memory, std::vector<symbolic_container> & symbolic_mapping, std::string & str) : current_arg_(0), memory_(memory), symbolic_mapping_(symbolic_mapping), str_(str) { }

          void call_on_leaf(std::size_t index, statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element) const {
            if(type_family==HOST_SCALAR_TYPE_FAMILY)
              header_value_generation(type_family,type,element, symbolic_mapping_[index]);
            else
              header_pointer_generation(type_family,type,element, symbolic_mapping_[index]);
          }

          void call_on_op(operation_node_type_family, operation_node_type type) const {  }

          void call_before_expansion() const { }

          void call_after_expansion() const { }
      };


    }

  }

}
#endif
