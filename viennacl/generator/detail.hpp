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

#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/result_of.hpp"
namespace viennacl{

  namespace generator{

    namespace detail{

      using namespace viennacl::scheduler;


      class mapped_container;

      template<class TraversalFunctor>
      void traverse(statement::container_type const & array, TraversalFunctor const & fun, bool use_extended_leaf = false, std::size_t index = 0);
      void generate(std::string const & index, utils::kernel_generation_stream & stream, mapped_container const & s);
      const char * generate(operation_node_type arg);
      const char * generate(statement_node_type arg);

      enum leaf_type{
        LHS_LEAF_TYPE,
        RHS_LEAF_TYPE
      };

      typedef std::pair<std::size_t, leaf_type> leaf_descriptor;
      typedef std::map< leaf_descriptor, tools::shared_ptr<detail::mapped_container> > mapping_type;

      struct expression_tree_ref{
          std::size_t node_index_;
          statement_node node_;
          statement::container_type const * array_;
          mapping_type const * mapping_;
      };

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
          std::string index_string_;
          utils::kernel_generation_stream & stream_;
          mapping_type const & mapping_;
        public:
          expression_generation_traversal(std::string const & index, utils::kernel_generation_stream & stream, mapping_type const & mapping) : index_string_(index), stream_(stream), mapping_(mapping){ }
          void call_on_leaf(std::size_t index, leaf_type lhs_rhs, statement_node const & node, statement::container_type const * array) const {
            generate(index_string_, stream_,*mapping_.at(std::make_pair(index, lhs_rhs)));
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

      /** @brief Base class for mapping viennacl datastructure to generator-friendly structures
       */
      class mapped_container{
        protected:
          std::string scalartype_;

        public:
          mapped_container(std::string const & scalartype) : scalartype_(scalartype){ }
          std::string const & scalartype() { return scalartype_; }
          virtual void fetch(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream){ }
          virtual void write_back(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream){ }
          virtual void generate(std::string const & index, utils::kernel_generation_stream & stream) const = 0;
          virtual ~mapped_container(){ }
      };

      /** @brief Mapping of a host scalar to a generator class */
      class mapped_host_scalar : public mapped_container{
          friend class prototype_generation_traversal;
          std::string name_;
        public:
          mapped_host_scalar(std::string const & scalartype) : mapped_container(scalartype){ }
          void generate(std::string const & index, utils::kernel_generation_stream & stream) const{
              stream << name_;
          }
      };

      /** @brief Base class for datastructures passed by pointer */
      class mapped_handle : public mapped_container{
        protected:
          std::string name_;
          std::string access_name_;

          virtual void offset(std::string const & index, utils::kernel_generation_stream & stream) const = 0;
        public:
          std::string const & name() { return name_; }
          void access_name(std::string const & str) { access_name_ = str; }
          std::string const & access_name() const { return access_name_; }
          mapped_handle(std::string const & scalartype) : mapped_container(scalartype){ }
          void generate(std::string const & index, utils::kernel_generation_stream & stream) const{
            if(!access_name_.empty())
              stream << access_name_;
            else{
              stream << name_;
              stream << '[' ;
              offset(index, stream);
              stream << ']';
            }
          }
      };

      /** @brief Mapping of a vector to a generator class */
      class mapped_vector : public mapped_handle{
          friend class prototype_generation_traversal;

          expression_tree_ref access_;

          std::string start_name_;
          std::string stride_name_;
          std::string shift_name_;

          void offset(std::string const & index, utils::kernel_generation_stream & stream) const {
            if(access_.mapping_==NULL)
              stream << index;
            else if(access_.node_.rhs_type_==COMPOSITE_OPERATION_TYPE)
              traverse(*access_.array_, expression_generation_traversal(index, stream,*access_.mapping_), false, access_.node_.rhs_.node_index_);
            else
              detail::generate(index, stream,*access_.mapping_->at(std::make_pair(access_.node_index_, RHS_LEAF_TYPE)));
          }

        public:
          mapped_vector(std::string const & scalartype) : mapped_handle(scalartype){ }
      };

      /** @brief Mapping of a matrix to a generator class */
      class mapped_matrix : public mapped_handle{
          friend class prototype_generation_traversal;

          expression_tree_ref access1_;
          std::string start1_name_;
          std::string stride1_name_;
          std::string shift1_name_;

          expression_tree_ref access2_;
          std::string start2_name_;
          std::string stride2_name_;
          std::string shift2_name_;

          bool is_row_major_;

          void offset(std::string const & index, utils::kernel_generation_stream & stream) const {
            stream << index;
          }

        public:
          bool is_row_major() const { return is_row_major_; }
          mapped_matrix(std::string const & scalartype) : mapped_handle(scalartype){ }
      };

      /** @brief Mapping of a symbolic vector to a generator class */
      class mapped_symbolic_vector : public mapped_container{
          friend class prototype_generation_traversal;
          std::string value_name_;
          std::string index_name_;
          bool is_value_static_;
        public:
          mapped_symbolic_vector(std::string const & scalartype) : mapped_container(scalartype){ }
          void generate(std::string const & index, utils::kernel_generation_stream & stream) const{
            stream << value_name_;
          }
      };

      /** @brief Mapping of a symbolic matrix to a generator class */
      class mapped_symbolic_matrix : public mapped_container{
          friend class prototype_generation_traversal;
          std::string value_name_;
          bool is_diag_;
        public:
          mapped_symbolic_matrix(std::string const & scalartype) : mapped_container(scalartype){ }
          void generate(std::string const & index, utils::kernel_generation_stream & stream) const{
            stream << value_name_;
          }
      };

      void generate(std::string const & index, utils::kernel_generation_stream & stream, mapped_container const & s){
        s.generate(index, stream);
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

      ///TODO : Add distinction among vector_base, symbolic_vector_base, ...
      const char * generate(statement_node_type type){
        switch(type){
          case COMPOSITE_OPERATION_TYPE : return "";

          //vector:
          case VECTOR_FLOAT_TYPE : return "vf";

          //symbolic vector:
          case SYMBOLIC_VECTOR_FLOAT_TYPE : return "svf";

         //matrix row:
         case MATRIX_ROW_FLOAT_TYPE : return "mrf";

         //matrix col:
         case MATRIX_COL_FLOAT_TYPE : return "mcf";

         //symbolic matrix:
         case SYMBOLIC_MATRIX_FLOAT_TYPE : return "smf";

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

          template<class ScalarType>
          void host_scalar_prototype(mapping_type::key_type const & key, ScalarType * scal) const {
            mapped_host_scalar * p = new mapped_host_scalar(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
            p->name_ = prototype_value_generation(p->scalartype_, (void *)scal);
          }

          template<class ScalarType>
          void vector_prototype(std::size_t index, leaf_type lhs_rhs, statement_node const & node, statement::container_type const * array, mapping_type::key_type const & key, vector_base<ScalarType> * vec) const {
            mapped_vector * p = new mapped_vector(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
            if(array){
              p->access_.node_index_ = index;
              p->access_.node_ = node;
              p->access_.array_ = array;
              p->access_.mapping_ = &mapping_;
            }
            p->name_ = prototype_pointer_generation(p->scalartype_, (void*)vec);

            if(vec->start() > 0)
              p->start_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
            if(vec->stride() > 1)
              p->shift_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
          }

          template<class ScalarType, class F>
          void matrix_prototype(mapping_type::key_type const & key, matrix_base<ScalarType, F> * mat) const {
            mapped_matrix * p = new mapped_matrix(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
            p->name_ = prototype_pointer_generation(p->scalartype_, (void*)mat);

            if(utils::is_same_type<F, viennacl::row_major>::value)
               p->is_row_major_ = true;
            else
              p->is_row_major_ = false;

            if(mat->start1() > 0)
              p->start1_name_ = prototype_value_generation(p->scalartype_, (void*)mat);
            if(mat->stride1() > 1)
              p->stride1_name_ = prototype_value_generation(p->scalartype_, (void*)mat);

            if(mat->start2() > 0)
              p->start2_name_ = prototype_value_generation(p->scalartype_, (void*)mat);
            if(mat->stride2() > 1)
              p->stride2_name_ = prototype_value_generation(p->scalartype_, (void*)mat);

          }


          template<class ScalarType>
          void symbolic_vector_prototype(mapping_type::key_type const & key, symbolic_vector_base<ScalarType> * vec) const {
            mapped_symbolic_vector * p = new mapped_symbolic_vector(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));

            if(!vec->is_value_static())
              p->value_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
            if(vec->index().first)
              p->index_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
          }

          template<class ScalarType>
          void symbolic_matrix_prototype(mapping_type::key_type const & key, symbolic_matrix_base<ScalarType> * mat) const {
            mapped_symbolic_matrix * p = new mapped_symbolic_matrix(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));

            if(!mat->is_value_static())
              p->value_name_ = prototype_value_generation(p->scalartype_, (void*)mat);
            if(mat->diag())
              p->is_diag_ = true;
          }

        public:

          prototype_generation_traversal(std::map<void*, std::size_t> & memory, mapping_type & mapping, std::string & str, std::size_t & current_arg) : memory_(memory), mapping_(mapping), str_(str), current_arg_(current_arg) { }

          void call_on_leaf(std::size_t index, leaf_type lhs_rhs, statement_node const & node, statement::container_type const * array) const {
            mapping_type::key_type key = std::make_pair(index, lhs_rhs);
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
            switch(type_family){
              case HOST_SCALAR_TYPE_FAMILY:
                switch(type){
                  case HOST_SCALAR_FLOAT_TYPE : host_scalar_prototype(key, &element.host_float_); break;
                  default : throw "not implemented";
                }
                break;
              case VECTOR_TYPE_FAMILY:
                switch(type){
                  case VECTOR_FLOAT_TYPE : vector_prototype(index, lhs_rhs, node, array, key, element.vector_float_); break;
                  default : throw "not implemented";
                }
                break;
              case SYMBOLIC_VECTOR_TYPE_FAMILY:
                switch(type){
                  case SYMBOLIC_VECTOR_FLOAT_TYPE : symbolic_vector_prototype(key, element.symbolic_vector_float_); break;
                  default : throw "not implemented";
                }
                break;
              case MATRIX_ROW_TYPE_FAMILY:
                switch(type){
                  case MATRIX_ROW_FLOAT_TYPE : matrix_prototype(key, element.matrix_row_float_); break;
                  default : throw "not implemented";
                }
                break;
              case MATRIX_COL_TYPE_FAMILY:
                switch(type){
                  case MATRIX_COL_FLOAT_TYPE : matrix_prototype(key, element.matrix_col_float_); break;
                  default : throw "not implemented";
                }
                break;
              case SYMBOLIC_MATRIX_TYPE_FAMILY:
                switch(type){
                  case SYMBOLIC_MATRIX_FLOAT_TYPE : symbolic_matrix_prototype(key, element.symbolic_matrix_float_); break;
                  default : throw "not implemented";
                }
                break;              default:
                throw "not implemented";

            }
          }
      };

    }

  }

}
#endif
