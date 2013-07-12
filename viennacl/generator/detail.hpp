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

#include "viennacl/vector.hpp"
#include "viennacl/tools/shared_ptr.hpp"

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

      enum leaf_type{
        LHS_LEAF_TYPE,
        RHS_LEAF_TYPE
      };

      typedef std::pair<std::size_t, leaf_type> leaf_descriptor;

      typedef std::map< leaf_descriptor, tools::shared_ptr<detail::symbolic_container> > mapping_type;

      std::string generate_value_kernel_argument(std::string const & scalartype, std::string const & name){
        return scalartype + ' ' + name + ",";
      }

      std::string generate_pointer_kernel_argument(std::string const & address_space, std::string const & scalartype, std::string const & name){
        return "__global " +  scalartype + "* " + name + ",";
      }

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
          void call_on_leaf(std::size_t index, leaf_type lhs_rhs, statement_node const & node, statement::container_type const * array) const {
            generate(stream_,*mapping_.at(std::make_pair(index, lhs_rhs)));
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

      struct host_scalar_descriptor{
          std::string name_;
      };

      class symbolic_container{
        private:
          std::string generate_scalartype(statement_node_type type){
            switch(type){
              case VECTOR_FLOAT_TYPE : return "float";

              case VECTOR_INITIALIZER_FLOAT_TYPE : return "float";
              default: throw "unrecognized type";
            }
          }
          virtual void generate_impl(utils::kernel_generation_stream & stream) const = 0;

        protected:
          std::string scalartype_;
          std::string access_name_;

        public:
          symbolic_container(statement_node_type type) : scalartype_(generate_scalartype(type)){ }
          virtual void fetch(std::set<std::string> & fetched, utils::kernel_generation_stream & stream){ }
          virtual void write_back(std::set<std::string> & fetched, utils::kernel_generation_stream & stream){ }
          void generate(utils::kernel_generation_stream & stream) const{
            if(!access_name_.empty())
              stream << access_name_;
            else
              generate_impl(stream);
          }
          virtual ~symbolic_container(){ }
      };



      class symbolic_host_scalar : public symbolic_container{
          friend class prototype_generation_traversal;
          std::string name_;

          void generate_impl(utils::kernel_generation_stream & stream) const{
              stream << name_;
          }
        public:
          symbolic_host_scalar(statement_node_type type) : symbolic_container(type){ }
      };

      template<class SCALARTYPE>
      class symbolic_vector : public symbolic_container{
          friend class prototype_generation_traversal;

          std::size_t node_index_;
          statement_node node_;
          statement::container_type const * array_;
          mapping_type const * mapping_;

          std::string name_;

          std::string start_name_;
          std::string stride_name_;
          std::string shift_name_;

          vector_base<SCALARTYPE> const * vec_;

          void generate_impl(utils::kernel_generation_stream & stream) const{
              stream << name_;
              stream << '[' ;
              if(mapping_==NULL)
                stream << "i";
              else if(node_.rhs_type_==COMPOSITE_OPERATION_TYPE)
                traverse(*array_, expression_generation_traversal(stream,*mapping_), false, node_.rhs_.node_index_);
              else
                detail::generate(stream,*mapping_->at(std::make_pair(node_index_, RHS_LEAF_TYPE)));
              stream << ']';
          }

        public:
          symbolic_vector(statement_node_type type, vector_base<SCALARTYPE> const * vec) : symbolic_container(type), vec_(vec){ }

          void fetch(std::set<std::string> & fetched, utils::kernel_generation_stream & stream){
            std::string new_access_name = name_ + "_private";
            if(fetched.find(name_)==fetched.end()){
              stream << scalartype_ << " " << new_access_name << " = ";
              generate(stream);
              stream << ';' << std::endl;
              fetched.insert(name_);
            }
            access_name_ = new_access_name;
          }

          void write_back(std::set<std::string> & fetched, utils::kernel_generation_stream & stream){
            std::string old_access_name = access_name_;
            access_name_.clear();
            if(fetched.find(name_)!=fetched.end()){
              generate(stream);
              stream << " = " << old_access_name << ';' << std::endl;
              fetched.erase(name_);
            }
          }
      };

      template<class SCALARTYPE>
      class symbolic_vector_initializer : public symbolic_container{
          friend class prototype_generation_traversal;
          std::string value_name_;
          std::string index_name_;
          vector_initializer_base<SCALARTYPE> * vec_;
          void generate_impl(utils::kernel_generation_stream & stream) const{
            stream << value_name_;
          }
        public:
          symbolic_vector_initializer(statement_node_type type, vector_initializer_base<SCALARTYPE> * vec) : symbolic_container(type), vec_(vec){ }
      };

      void generate(utils::kernel_generation_stream & stream, symbolic_container const & s){
        s.generate(stream);
      }


      const char * generate(operation_node_type type){
        // unary expression
        switch(type){
          case OPERATION_UNARY_ABS_TYPE : return "abs";

          case OPERATION_BINARY_ASSIGN_TYPE : return "=";
          case OPERATION_BINARY_ADD_TYPE : return "+";
          case OPERATION_BINARY_ACCESS : return "";

          default : throw "not implemented";
        }
      }

      const char * generate(statement_node_type type){
        switch(type){
          case COMPOSITE_OPERATION_TYPE : return "";

          case VECTOR_FLOAT_TYPE : return "vf";
          default : throw "not implemented";
        }
      }


      template<class TraversalFunctor>
      void traverse(statement::container_type const & array, TraversalFunctor const & fun, bool stop_recursion_on_access = false, std::size_t index = 0){
        statement::value_type const & element = array[index];
        if(element.op_family_==OPERATION_UNARY_TYPE_FAMILY){
          fun.call_on_op(element.op_family_, element.op_type_);
          fun.call_before_expansion();
          if(element.lhs_type_==COMPOSITE_OPERATION_TYPE)
            traverse(array, fun, stop_recursion_on_access, element.lhs_.node_index_);
          else
            fun.call_on_leaf(index, LHS_LEAF_TYPE, element, NULL);
          fun.call_after_expansion();
        }
        if(element.op_family_==OPERATION_BINARY_TYPE_FAMILY){
          fun.call_before_expansion();
            if(element.lhs_type_==COMPOSITE_OPERATION_TYPE)
              traverse(array, fun, stop_recursion_on_access, element.lhs_.node_index_);
            else{
              if(element.op_type_==OPERATION_BINARY_ACCESS)
                fun.call_on_leaf(index, LHS_LEAF_TYPE, element, &array);
              else
                fun.call_on_leaf(index, LHS_LEAF_TYPE, element, NULL);
            }
            if(stop_recursion_on_access && element.op_type_==OPERATION_BINARY_ACCESS)
              return;

            fun.call_on_op(element.op_family_, element.op_type_);
            if(element.rhs_type_==COMPOSITE_OPERATION_TYPE)
              traverse(array, fun, stop_recursion_on_access, element.rhs_.node_index_);
            else
              fun.call_on_leaf(index, RHS_LEAF_TYPE, element, NULL);
            fun.call_after_expansion();
        }
      }




      class name_generation_traversal : public traversal_functor{
          std::string & str_;
        public:
          name_generation_traversal(std::string & str) : str_(str) { }
          void call_on_leaf(std::size_t index, leaf_type lhs_rhs, statement_node const & node, statement::container_type const * array) const {
            if(lhs_rhs == LHS_LEAF_TYPE)
              str_ += detail::generate(node.lhs_type_);
            else
              str_ += detail::generate(node.rhs_type_);
          }
          void call_on_op(operation_node_type_family, operation_node_type type) const { str_ += detail::generate(type); }
          void call_before_expansion() const { str_ += '('; }
          void call_after_expansion() const { str_ += ')'; }
      };

      class prototype_generation_traversal : public traversal_functor{
          std::map<void *, std::size_t> & memory_;
          mapping_type & mapping_;
          std::string & str_;
          std::size_t & current_arg_;


          std::string prototype_value_generation(std::string const & scalartype, void * handle) const{
            if(memory_.insert(std::make_pair(handle, current_arg_)).second){
              std::string name =  "arg" + utils::to_string(current_arg_++);
              str_ += detail::generate_value_kernel_argument(scalartype, name);
              return name;
            }
            else
              return "arg" + utils::to_string(memory_.at(handle));
          }

          std::string prototype_pointer_generation(std::string const & scalartype, void * handle) const {
            if(memory_.insert(std::make_pair(handle, current_arg_)).second){
              std::string name =  "arg" + utils::to_string(current_arg_++);
              str_ += detail::generate_pointer_kernel_argument("__global", scalartype, name);
              return name;
            }
            else
              return "arg" + utils::to_string(memory_.at(handle));
          }

          void host_scalar_prototype(symbolic_host_scalar * p) const {
            p->name_ = prototype_value_generation(p->scalartype_, (void *)p);
          }

          template<typename T>
          void vector_initializer_prototype(symbolic_vector_initializer<T> * p) const {
            if(!p->vec_->is_value_static())
              p->value_name_ = prototype_value_generation(p->scalartype_, (void *)p->vec_);
            if(p->vec_->index().first)
              p->index_name_ = prototype_value_generation(p->scalartype_, (void *)p->vec_);
          }

          template<typename T>
          void vector_prototype(symbolic_vector<T> * p) const {
            p->name_ = prototype_pointer_generation(p->scalartype_, (void *)p->vec_);
            if(p->vec_->start() > 0)
              p->start_name_ = prototype_value_generation(p->scalartype_, (void *)p->vec_);
            if(p->vec_->stride() > 1)
              p->shift_name_ = prototype_value_generation(p->scalartype_, (void *)p->vec_);
          }

        public:

          prototype_generation_traversal(std::map<void*, std::size_t> & memory, mapping_type & mapping, std::string & str, std::size_t & current_arg) : memory_(memory), mapping_(mapping), str_(str), current_arg_(current_arg) { }

          void call_on_leaf(std::size_t index, leaf_type lhs_rhs, statement_node const & node, statement::container_type const * array) const {
            statement_node_type type;
            statement_node_type_family type_family;
            lhs_rhs_element element;
            if(lhs_rhs==LHS_LEAF_TYPE){
               type = node.lhs_type_;
               type_family = node.lhs_type_family_;
               element = node.lhs_;
            }
            else{
              type = node.rhs_type_;
              type_family = node.rhs_type_family_;
              element = node.rhs_;
            }
            if(type_family==HOST_SCALAR_TYPE_FAMILY){
              symbolic_host_scalar * p = new symbolic_host_scalar(type);
              mapping_.insert(std::make_pair(std::make_pair(index, lhs_rhs), tools::shared_ptr<symbolic_container>(p)));
              host_scalar_prototype(p);
            }
            else if(type_family==VECTOR_TYPE_FAMILY){
              if(type==VECTOR_FLOAT_TYPE){
                symbolic_vector<float> * p = new symbolic_vector<float>(type, element.vector_float_);
                mapping_.insert(std::make_pair(std::make_pair(index, lhs_rhs), tools::shared_ptr<symbolic_container>(p)));
                if(array){
                  p->node_index_ = index;
                  p->node_ = node;
                  p->array_ = array;
                  p->mapping_ = &mapping_;
                }
                vector_prototype(p);
              }
              if(type==VECTOR_DOUBLE_TYPE){
                symbolic_vector<double> * p = new symbolic_vector<double>(type, element.vector_double_);
                mapping_.insert(std::make_pair(std::make_pair(index, lhs_rhs), tools::shared_ptr<symbolic_container>(p)));
                if(array){
                  p->node_index_ = index;
                  p->node_ = node;
                  p->array_ = array;
                  p->mapping_ = &mapping_;
                }
                vector_prototype(p);
              }
            }
            else if(type_family==VECTOR_INITIALIZER_TYPE_FAMILY){
              if(type==VECTOR_INITIALIZER_FLOAT_TYPE){
                symbolic_vector_initializer<float> * p = new symbolic_vector_initializer<float>(type, element.vector_initializer_float_);
                mapping_.insert(std::make_pair(std::make_pair(index, lhs_rhs), tools::shared_ptr<symbolic_container>(p)));
                vector_initializer_prototype(p);
              }
              if(type==VECTOR_INITIALIZER_DOUBLE_TYPE){
                symbolic_vector_initializer<double> * p = new symbolic_vector_initializer<double>(type, element.vector_initializer_double_);
                mapping_.insert(std::make_pair(std::make_pair(index, lhs_rhs), tools::shared_ptr<symbolic_container>(p)));
                vector_initializer_prototype(p);
              }
            }
          }
      };

    }

  }

}
#endif
