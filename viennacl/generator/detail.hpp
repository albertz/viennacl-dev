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

      enum node_type{
        LHS_NODE_TYPE,
        PARENT_TYPE,
        RHS_NODE_TYPE
      };
      typedef std::pair<std::size_t, node_type> index_info;
      class mapped_container;
      typedef std::map< index_info, tools::shared_ptr<detail::mapped_container> > mapping_type;


      template<class TraversalFunctor>
      void traverse(statement::container_type const & array, TraversalFunctor const & fun, bool shallow_traversal = false, index_info const & key = std::make_pair(0, PARENT_TYPE));
      std::string generate(std::string const & index, mapped_container const & s);
      const char * generate(operation_node_type arg);
      const char * generate(statement_node_type arg);


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
          void call_on_leaf(index_info const & key, statement_node const & node,  statement::container_type const * array) const {
            stream_ << generate(index_string_,*mapping_.at(key));
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


      /** @brief Base class for mapping viennacl datastructure to generator-friendly structures
       */
      class mapped_container{
        protected:
          struct node_info{
              statement::container_type const * array_;
              index_info index_;
          };

        protected:
          std::string scalartype_;

        public:
          mapped_container(std::string const & scalartype) : scalartype_(scalartype){ }
          std::string const & scalartype() { return scalartype_; }
          virtual void fetch(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream){ }
          virtual void write_back(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream){ }
          virtual std::string generate(std::string const & index) const = 0;
          virtual ~mapped_container(){ }
      };

      /** @brief Mapping of a host scalar to a generator class */
      class mapped_host_scalar : public mapped_container{
          friend class map_generate_prototype;
          std::string name_;
        public:
          mapped_host_scalar(std::string const & scalartype) : mapped_container(scalartype){ }
          std::string generate(std::string const & index) const{
              return name_;
          }
      };

      /** @brief Base class for datastructures passed by pointer */
      class mapped_handle : public mapped_container{
        protected:
          std::string name_;
          std::string access_name_;

          virtual std::string offset(std::string const & index) const = 0;
        public:
          std::string const & name() { return name_; }
          void access_name(std::string const & str) { access_name_ = str; }
          std::string const & access_name() const { return access_name_; }
          mapped_handle(std::string const & scalartype) : mapped_container(scalartype){ }

          void fetch(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) {
            std::string new_access_name = name_ + "_private";
            if(fetched.find(name_)==fetched.end()){
              stream << scalartype_ << " " << new_access_name << " = " << generate(index) << ';' << std::endl;
              fetched.insert(name_);
            }
            access_name_ = new_access_name;
          }

          void write_back(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) {
            std::string old_access_name = access_name_ ;
            access_name_ = "";
            if(fetched.find(name_)!=fetched.end()){
              stream << generate(index) << " = " << old_access_name << ';' << std::endl;
              fetched.erase(name_);
            }
          }

          std::string generate(std::string const & index) const{
            if(!access_name_.empty())
              return access_name_;
            else
              return name_  + '[' + offset(index) + ']';
          }
      };


      /** @brief Mapping of a vector to a generator class */
      class mapped_vector : public mapped_handle{
          friend class map_generate_prototype;

          node_info access_node_;

          std::string start_name_;
          std::string stride_name_;
          std::string shift_name_;

          std::string offset(std::string const & index) const {
              return index;
          }

        public:
          mapped_vector(std::string const & scalartype) : mapped_handle(scalartype){ }
      };

      /** @brief Mapping of a matrix to a generator class */
      class mapped_matrix : public mapped_handle{
          friend class map_generate_prototype;

          std::string start1_name_;
          std::string stride1_name_;
          std::string shift1_name_;

          std::string start2_name_;
          std::string stride2_name_;
          std::string shift2_name_;

          bool is_row_major_;

          std::string offset(std::string const & index) const {
            return index;
          }

        public:
          bool is_row_major() const { return is_row_major_; }
          mapped_matrix(std::string const & scalartype) : mapped_handle(scalartype){ }
      };

      /** @brief Mapping of a symbolic vector to a generator class */
      class mapped_symbolic_vector : public mapped_container{
          friend class map_generate_prototype;
          std::string value_name_;
          std::string index_name_;
          bool is_value_static_;
        public:
          mapped_symbolic_vector(std::string const & scalartype) : mapped_container(scalartype){ }
          std::string generate(std::string const & index) const{
            return value_name_;
          }
      };

      /** @brief Mapping of a symbolic matrix to a generator class */
      class mapped_symbolic_matrix : public mapped_container{
          friend class map_generate_prototype;
          std::string value_name_;
          bool is_diag_;
        public:
          mapped_symbolic_matrix(std::string const & scalartype) : mapped_container(scalartype){ }
          std::string generate(std::string const & index) const{
            return value_name_;
          }
      };

      std::string generate(std::string const & index, mapped_container const & s){
        return s.generate(index);
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

      index_info get_new_key(statement_node_type_family type_family, std::size_t current_index, std::size_t next_index, node_type node_tag){
        if(type_family==COMPOSITE_OPERATION_FAMILY)
          return std::make_pair(next_index, PARENT_TYPE);
        else
          return std::make_pair(current_index, node_tag);
      }


      template<class TraversalFunctor>
      void traverse(statement::container_type const & array, TraversalFunctor const & fun, bool shallow_traversal, index_info const & key){
        std::size_t index = key.first;
        std::size_t node_tag = key.second;
        statement::value_type const & element = array[index];
        if(node_tag == PARENT_TYPE){
          if(element.op_family_==OPERATION_UNARY_TYPE_FAMILY){
            fun.call_on_op(element.op_family_, element.op_type_);
            fun.call_before_expansion();
            traverse(array, fun, shallow_traversal, get_new_key(element.lhs_type_family_, index, element.lhs_.node_index_, LHS_NODE_TYPE));
            fun.call_after_expansion();
          }
          if(element.op_family_==OPERATION_BINARY_TYPE_FAMILY){
            fun.call_before_expansion();
            traverse(array, fun, shallow_traversal, get_new_key(element.lhs_type_family_, index, element.lhs_.node_index_, LHS_NODE_TYPE));
            fun.call_on_op(element.op_family_, element.op_type_);
            traverse(array, fun, shallow_traversal, get_new_key(element.rhs_type_family_, index, element.rhs_.node_index_, RHS_NODE_TYPE));
            fun.call_after_expansion();
          }
        }
        else{
          fun.call_on_leaf(key, element, &array);
        }
      }

//      class name_generation_traversal : public traversal_functor{
//          std::string & str_;
//        public:
//          name_generation_traversal(std::string & str) : str_(str) { }
//          void call_on_leaf(std::size_t index, node_type lhs_rhs, statement_node const & node) const {
//            if(lhs_rhs == LHS_NODE_TYPE)
//              str_ += detail::generate(node.lhs_type_);
//            else
//              str_ += detail::generate(node.rhs_type_);
//          }
//          void call_on_op(operation_node_type_family, operation_node_type type) const { str_ += detail::generate(type); }
//          void call_before_expansion() const { str_ += '('; }
//          void call_after_expansion() const { str_ += ')'; }
//      };

      class map_generate_prototype : public traversal_functor{
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
          mapped_host_scalar * host_scalar_prototype(index_info const & key, ScalarType * scal) const {
            mapped_host_scalar * p = new mapped_host_scalar(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
            p->name_ = prototype_value_generation(p->scalartype_, (void *)scal);
            return p;
          }

          template<class ScalarType>
          mapped_vector * vector_prototype(index_info const & key, vector_base<ScalarType> * vec) const {
            mapped_vector * p = new mapped_vector(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
            p->name_ = prototype_pointer_generation(p->scalartype_, (void*)vec);
            if(vec->start() > 0)
              p->start_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
            if(vec->stride() > 1)
              p->shift_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
            return p;
          }

          template<class ScalarType, class F>
          mapped_matrix * matrix_prototype(index_info const & key, matrix_base<ScalarType, F> * mat) const {
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
            return p;
          }


          template<class ScalarType>
          mapped_symbolic_vector * symbolic_vector_prototype(index_info const & key, symbolic_vector_base<ScalarType> * vec) const {
            mapped_symbolic_vector * p = new mapped_symbolic_vector(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
            if(!vec->is_value_static())
              p->value_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
            if(vec->index().first)
              p->index_name_ = prototype_value_generation(p->scalartype_, (void*)vec);
            return p;
          }

          template<class ScalarType>
          mapped_symbolic_matrix * symbolic_matrix_prototype(index_info const & key, symbolic_matrix_base<ScalarType> * mat) const {
            mapped_symbolic_matrix * p = new mapped_symbolic_matrix(utils::type_to_string<ScalarType>::value());
            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
            if(!mat->is_value_static())
              p->value_name_ = prototype_value_generation(p->scalartype_, (void*)mat);
            if(mat->diag())
              p->is_diag_ = true;
            return p;
          }

        public:

          map_generate_prototype(std::map<void*, std::size_t> & memory, mapping_type & mapping, std::string & str, std::size_t & current_arg) : memory_(memory), mapping_(mapping), str_(str), current_arg_(current_arg) { }

          void call_on_leaf(index_info const & key, statement_node const & node,  statement::container_type const * array) const {
            statement_node_type type;
            statement_node_type_family type_family;
            lhs_rhs_element element;
            if(key.second==LHS_NODE_TYPE){
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
              {
                switch(type){
                  case HOST_SCALAR_FLOAT_TYPE : host_scalar_prototype(key, &element.host_float_); break;
                  default : throw "not implemented";
                }
                break;
              }
              case VECTOR_TYPE_FAMILY:
              {
                mapped_vector * v = NULL;
                switch(type){
                  case VECTOR_FLOAT_TYPE : v = vector_prototype(key, element.vector_float_);  break;
                  default : throw "not implemented";
                }
                if(node.op_type_==OPERATION_BINARY_ACCESS){
                  v->access_node_.array_ = array;
                  v->access_node_.index_ = get_new_key(node.rhs_type_family_, key.first, node.rhs_.node_index_, RHS_NODE_TYPE);
                }
                break;
              }
              case SYMBOLIC_VECTOR_TYPE_FAMILY:
              {
                switch(type){
                  case SYMBOLIC_VECTOR_FLOAT_TYPE : symbolic_vector_prototype(key, element.symbolic_vector_float_); break;
                  default : throw "not implemented";
                }
                break;
              }
              case MATRIX_ROW_TYPE_FAMILY:
              {
                switch(type){
                  case MATRIX_ROW_FLOAT_TYPE : matrix_prototype(key, element.matrix_row_float_); break;
                  default : throw "not implemented";
                }
                break;
              }
              case MATRIX_COL_TYPE_FAMILY:
              {
                switch(type){
                  case MATRIX_COL_FLOAT_TYPE : matrix_prototype(key, element.matrix_col_float_); break;
                  default : throw "not implemented";
                }
                break;
              }
              case SYMBOLIC_MATRIX_TYPE_FAMILY:
              {
                switch(type){
                  case SYMBOLIC_MATRIX_FLOAT_TYPE : symbolic_matrix_prototype(key, element.symbolic_matrix_float_); break;
                  default : throw "not implemented";
                }
                break;
              }
              default : throw "not implemented";
            }
          }

      };

    }

  }

}
#endif
