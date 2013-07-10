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

#include <set>

#include "CL/cl.h"

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/utils.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{

      using namespace viennacl::scheduler;

      class symbolic_container;

      template<class TraversalFunctor>
      void traverse(statement::container_type const & array, TraversalFunctor const & fun, bool use_extended_leaf = false, std::size_t index = 0);
      void generate(utils::kernel_generation_stream & stream, symbolic_container const & s);
      const char * generate(operation_node_type arg);
      const char * generate(statement_node_type arg);

      typedef std::map<std::size_t, detail::symbolic_container> mapping_type;

      class traversal_functor{
        public:
          void call_on_op(operation_node_type_family, operation_node_type type) const { }
          void call_before_expansion() const { }
          void call_after_expansion() const { }
      };

      class expression_generation_traversal : public traversal_functor{

        private:
          utils::kernel_generation_stream & stream_;
          mapping_type const & mapping_;

        public:
          expression_generation_traversal(utils::kernel_generation_stream & stream, mapping_type const & mapping) : stream_(stream), mapping_(mapping){ }
          void call_on_extended_leaf(std::size_t index, statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element, statement::container_type const &, std::size_t) const {
            call_on_leaf(index,type_family,type,element);
          }
          void call_on_leaf(std::size_t index, statement_node_type_family, statement_node_type, lhs_rhs_element) const {
            generate(stream_,mapping_.at(index));
          }
          void call_on_op(operation_node_type_family, operation_node_type type) const {
            stream_ << detail::generate(type);
          }
          void call_before_expansion() const {
            stream_ << '(';
          }
          void call_after_expansion() const {
            stream_ << ')';
          }
      };

      class symbolic_container{
          friend class prototype_generation_traversal;

        private:
          statement_node_type_family type_family_;
          statement_node_type type_;

          std::size_t access_index_;
          statement::container_type const * array_;
          mapping_type const * mapping_;

          lhs_rhs_element element_;
          std::string scalartype_;
          std::string name_;
          std::string access_name_;

        private:
          void offset(utils::kernel_generation_stream & stream) const {
            if(array_==NULL)
              stream << "i";
            else
              traverse(*array_, expression_generation_traversal(stream,*mapping_), false, access_index_);
          }

        public:
          void fetch(std::set<std::string> & fetched, utils::kernel_generation_stream & stream){
            if(fetched.find(name_)==fetched.end()){
              std::string new_access_name = name_ + "_private";
              stream << scalartype_ << " " << new_access_name << " = ";
              generate(stream);
              stream << ';' << std::endl;
              access_name_ = new_access_name;
              fetched.insert(name_);
            }
          }

          void write_back(std::set<std::string> & fetched, utils::kernel_generation_stream & stream){
            if(fetched.find(name_)!=fetched.end()){
              std::string old_access_name = access_name_;
              access_name_.clear();
              generate(stream);
              stream << " = " << old_access_name << ';' << std::endl;
              fetched.erase(name_);
            }
          }

          void generate(utils::kernel_generation_stream & stream) const{
            if(!access_name_.empty())
              stream << access_name_;
            if(type_family_==HOST_SCALAR_TYPE_FAMILY)
              stream << name_;
            if(type_family_==SCALAR_TYPE_FAMILY)
              stream << '*' << name_;
            stream << name_ << "[" ;
            offset(stream);
            stream << "]";
          }
      };


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

      void generate(utils::kernel_generation_stream & stream, symbolic_container const & s){
        s.generate(stream);
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
      void traverse(statement::container_type const & array, TraversalFunctor const & fun, bool use_extended_leaf = false, std::size_t index = 0){
        statement::value_type const & element = array[index];
        if(element.op_family_==OPERATION_UNARY_TYPE_FAMILY){
          fun.call_on_op(element.op_family_, element.op_type_);
          fun.call_before_expansion();
          if(element.lhs_type_==COMPOSITE_OPERATION_TYPE)
            traverse(array, fun, use_extended_leaf, element.lhs_.node_index_);
          else
            fun.call_on_leaf(index, element.lhs_type_family_, element.lhs_type_, element.lhs_);
          fun.call_after_expansion();
        }
        if(element.op_family_==OPERATION_BINARY_TYPE_FAMILY){
          fun.call_before_expansion();
          if(use_extended_leaf && element.op_type_==OPERATION_BINARY_ACCESS){
            fun.call_on_extended_leaf(index, element.lhs_type_family_, element.lhs_type_, element.lhs_, array, element.rhs_.node_index_);
          }
          else{
            if(element.lhs_type_==COMPOSITE_OPERATION_TYPE)
              traverse(array, fun, use_extended_leaf, element.lhs_.node_index_);
            else
              fun.call_on_leaf(index, element.lhs_type_family_, element.lhs_type_, element.lhs_);
            fun.call_on_op(element.op_family_, element.op_type_);
            if(element.rhs_type_==COMPOSITE_OPERATION_TYPE)
              traverse(array, fun, use_extended_leaf, element.rhs_.node_index_);
            else
              fun.call_on_leaf(index, element.rhs_type_family_, element.rhs_type_, element.rhs_);
            fun.call_after_expansion();
          }
        }
      }




      class name_generation_traversal : public traversal_functor{
          std::string & str_;
        public:
          name_generation_traversal(std::string & str) : str_(str) { }
          void call_on_extended_leaf(std::size_t index, statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element, statement::container_type const &, std::size_t) const { call_on_leaf(index,type_family,type,element); }
          void call_on_leaf(std::size_t, statement_node_type_family, statement_node_type type, lhs_rhs_element) const { str_ += detail::generate(type); }
          void call_on_op(operation_node_type_family, operation_node_type type) const { str_ += detail::generate(type); }
          void call_before_expansion() const { str_ += '('; }
          void call_after_expansion() const { str_ += ')'; }
      };

      class prototype_generation_traversal : public traversal_functor{

          mutable unsigned int current_arg_;
          std::map<cl_mem, std::size_t> & memory_;
          mapping_type & mapping_;
          std::string & str_;


          void prototype_value_generation(statement_node_type_family, statement_node_type type, lhs_rhs_element, symbolic_container & sym) const{
            sym.name_ = "arg" + utils::to_string(current_arg_++);
            str_ += sym.scalartype_ + ' '  + sym.name_ + ",";
          }

          void prototype_pointer_generation(statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element, symbolic_container & sym) const {
            if(memory_.insert(std::make_pair(detail::get_handle(type, element), current_arg_)).second){
              sym.name_ =  "arg" + utils::to_string(current_arg_++);
              str_ += "__global " +  sym.scalartype_ + "* " + sym.name_ + ",";
            }
            else
              sym.name_ = "arg" + utils::to_string(memory_.at(detail::get_handle(type, element)));
          }

        public:
          prototype_generation_traversal(std::map<cl_mem, std::size_t> & memory, mapping_type & mapping, std::string & str) : current_arg_(0), memory_(memory), mapping_(mapping), str_(str) { }

          void call_on_extended_leaf(std::size_t index, statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element, statement::container_type const & array, std::size_t access_index) const {
            mapping_[index].array_ = &array;
            mapping_[index].mapping_ = &mapping_;
            mapping_[index].access_index_ = access_index;
            call_on_leaf(index, type_family, type, element);
          }

          void call_on_leaf(std::size_t index, statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element) const {
            mapping_[index].scalartype_ = detail::generate_scalartype(type);
            mapping_[index].type_ = type;
            mapping_[index].type_family_ = type_family;
            mapping_[index].element_ = element;
            if(type_family==HOST_SCALAR_TYPE_FAMILY)
              prototype_value_generation(type_family,type,element, mapping_[index]);
            else
              prototype_pointer_generation(type_family,type,element, mapping_[index]);
          }
      };

    }

  }

}
#endif
