#ifndef VIENNACL_GENERATOR_MAP_GENERATE_PROTOTYPE_HPP
#define VIENNACL_GENERATOR_MAP_GENERATE_PROTOTYPE_HPP

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


/** @file viennacl/generator/map_generate_prototype.hpp
    @brief Functor to map a statement and generate the prototype
*/

#include <set>

#include "viennacl/forwards.h"
#include "viennacl/scheduler/forwards.h"
#include "viennacl/generator/forwards.h"

#include "viennacl/tools/shared_ptr.hpp"

#include "viennacl/generator/generate_utils.hpp"
#include "viennacl/generator/utils.hpp"
#include "viennacl/generator/mapped_types.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{
      class map_functor{
          std::string create_name(unsigned int & current_arg, std::map<void *, std::size_t> & memory, void * handle) const{
            if(memory.insert(std::make_pair(handle, current_arg)).second)
              return "arg" + utils::to_string(current_arg_++);
            else
              return "arg" + utils::to_string(memory.at(handle));
          }

        public:
          map_functor(std::map<void *, std::size_t> & memory, unsigned int & current_arg) : memory_(memory), current_arg_(current_arg){ }

          //Binary leaf
          template<class T>
          result_type binary_leaf(unsigned int i, statement_node const & node,  statement::container_type const * array){
            T * p = new T("float");
            p->lhs_.array_ = &expr;
            p->lhs_.index_ = get_new_key(node.lhs_type_family_, i, node.lhs_.node_index_, LHS_NODE_TYPE);
            p->lhs_.mapping_ = &mapping_;

            p->rhs_.array_ = array;
            p->rhs_.index_ = get_new_key(node.rhs_type_family_, i, node.rhs_.node_index_, RHS_NODE_TYPE);
            p->rhs_.mapping_ = &mapping_;

            return p;
          }

          //Scalar mapping
          template<class T>
          result_type operator()(scalar<ScalarType> * scal){
            mapped_scalar * p = new mapped_scalar(utils::type_to_string<ScalarType>::value());
            p->name_ = create_name(current_arg_, memory_, (void*)scal);
          }

          //Vector mapping
          template<class T>
          result_type operator()(vector_base<ScalarType> * vec){
            mapped_vector * p = new mapped_vector(utils::type_to_string<ScalarType>::value());
            p->name_ = create_name(current_arg_, memory_, (void*)vec);
            if(vec->start() > 0)
              p->start_name_ = p->name_ +"_start";
            if(vec->stride() > 1)
              p->stride_name_ = p->name_ + "_stride";
            return p;
          }

        private:
          std::map<void *, std::size_t> & memory_;
          unsigned int & current_arg_;
          typedef mapped_vector * result_type;
      };

      void map_statement(scheduler::statement const & statement, std::map<void *, std::size_t> & memory, std::size_t & current_arg, mapping_type & mapping){
          scheduler::statement::container_type expr = statement.array();
          for(std::size_t i = 0 ; i < expr.size() ; ++i){
            scheduler::statement_node node = expr[i];

            if(node.lhs_type_family_!=COMPOSITE_OPERATION_FAMILY)
              mapping.insert(std::make_pair(std::make_pair(i, LHS_NODE_TYPE), utils::call_on_element(node.lhs_type_family_, node.lhs_type_, node.lhs_, map_functor(memory, current_arg))));

            if(node.op_type_==OPERATION_BINARY_INNER_PROD_TYPE){
              mapping.insert(std::make_pair(std::make_pair(i, PARENT_TYPE), map_functor(memory, current_arg).template binary_leaf<mapped_scalar_reduction>(i, node, &expr)));
            }

            if(node.rhs_type_family_!=COMPOSITE_OPERATION_FAMILY)
              mapping.insert(std::make_pair(std::make_pair(i, RHS_NODE_TYPE), utils::call_on_element(node.rhs_type_family_, node.rhs_type_, node.rhs_, map_functor(memory, current_arg))));
          }
      }

      template<class InputIterator>
      void map_all_statements(InputIterator begin, InputIterator end, mapping_type & mapping){
        std::map<void *, std::size_t> memory;
        std::size_t current_arg = 0;
        while(begin!=end){
          map_statement(*it, memory, current_arg, mapping);
        }
      }

//      class map_generate_prototype : public traversal_functor{
//          std::map<void *, std::size_t> & memory_;
//          mapping_type & mapping_;
//          std::string & str_;
//          std::size_t & current_arg_;


//          std::string create_name(void * handle) const{
//            if(memory_.insert(std::make_pair(handle, current_arg_)).second)
//              return "arg" + utils::to_string(current_arg_++);
//            else
//              return "arg" + utils::to_string(memory_.at(handle));
//          }

////          std::string create_name(void * handle) const {
////            if(memory_.insert(std::make_pair(handle, current_arg_)).second){
////              std::string name =  "arg" + utils::to_string(current_arg_++);
////              ;
////              return name;
////            }
////            else
////              return "arg" + utils::to_string(memory_.at(handle));
////          }

//          template<class ScalarType>
//          void host_scalar_prototype(index_info const & key, ScalarType * scal) const {
//            mapped_host_scalar * p = new mapped_host_scalar(utils::type_to_string<ScalarType>::value());
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            p->name_ = create_name((void *)scal);
//            str_ += detail::generate_value_kernel_argument(p->scalartype_, p->name_);
//          }

//          template<class ScalarType>
//          void scalar_prototype(index_info const & key, scalar<ScalarType> * scal) const{
//            mapped_scalar * p = new mapped_scalar(utils::type_to_string<ScalarType>::value());
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            p->name_ = create_name((void*)scal);
//            str_ += detail::generate_value_kernel_argument(p->scalartype_, p->name_);
//          }

//          template<class ScalarType>
//          void vector_prototype(index_info const & key, vector_base<ScalarType> * vec) const {
//            mapped_vector * p = new mapped_vector(utils::type_to_string<ScalarType>::value());
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            p->name_ = create_name((void*)vec);
//            str_ += detail::generate_pointer_kernel_argument("__global",p->scalartype_, p->name_);
//            if(vec->start() > 0){
//              p->start_name_ = p->name_ +"_start";
//              str_ += detail::generate_value_kernel_argument("unsigned int", p->start_name_);
//            }
//            if(vec->stride() > 1){
//              p->stride_name_ = p->name_ + "_stride";
//              str_ += detail::generate_value_kernel_argument("unsigned int", p->stride_name_);
//            }
//          }

//          template<class ScalarType, class F>
//          void matrix_prototype(index_info const & key, matrix_base<ScalarType, F> * mat) const {
//            mapped_matrix * p = new mapped_matrix(utils::type_to_string<ScalarType>::value());
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            p->name_ = create_name((void*)mat);
//            str_ += detail::generate_pointer_kernel_argument("__global",p->scalartype_, p->name_);
//            if(utils::is_same_type<F, viennacl::row_major>::value)
//               p->is_row_major_ = true;
//            else
//              p->is_row_major_ = false;
//            if(mat->start1() > 0){
//              p->start1_name_ = p->name_ + "start1_";
//              str_ += detail::generate_value_kernel_argument("unsigned int", p->start1_name_);
//            }
//            if(mat->stride1() > 1){
//              p->stride1_name_ = p->name_ + "stride1_";
//              str_ += detail::generate_value_kernel_argument("unsigned int", p->stride1_name_);
//            }
//            if(mat->start2() > 0){
//              p->start2_name_ = p->name_ + "start2_";
//              str_ += detail::generate_value_kernel_argument("unsigned int", p->start2_name_);
//            }
//            if(mat->stride2() > 1){
//              p->stride2_name_ = p->name_ + "stride2_";
//              str_ += detail::generate_value_kernel_argument("unsigned int", p->stride2_name_);
//            }

//          }

//          template<class ScalarType>
//          void symbolic_vector_prototype(index_info const & key, symbolic_vector_base<ScalarType> * vec) const {
//            mapped_symbolic_vector * p = new mapped_symbolic_vector(utils::type_to_string<ScalarType>::value());
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            if(!vec->is_value_static()){
//              p->value_name_ = create_name((void*)vec);
//              str_ += detail::generate_value_kernel_argument(p->scalartype(), p->value_name_);
//            }
//            if(vec->index().first)
//              p->index_name_ = create_name((void*)vec);
//            str_ += detail::generate_value_kernel_argument("unsigned int", p->index_name_);
//          }

//          template<class ScalarType>
//          void symbolic_matrix_prototype(index_info const & key, symbolic_matrix_base<ScalarType> * mat) const {
//            mapped_symbolic_matrix * p = new mapped_symbolic_matrix(utils::type_to_string<ScalarType>::value());
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            if(!mat->is_value_static()){
//              p->value_name_ = create_name((void*)mat);
//              str_ += detail::generate_value_kernel_argument(p->scalartype(), p->value_name_);
//            }
//            if(mat->diag())
//              p->is_diag_ = true;
//          }


//          void vector_reduction_prototype(index_info const & key, statement_node const & node,  statement::container_type const * array) const{
//            mapped_vector_reduction * p = new mapped_vector_reduction("float");
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            p->lhs_.array_ = array;
//            p->lhs_.index_ = get_new_key(node.lhs_type_family_, key.first, node.lhs_.node_index_, LHS_NODE_TYPE);
//            p->lhs_.mapping_ = &mapping_;
////            p->is_lhs_transposed_ = array->at(node.lhs_.node_index_).op_type_ == scheduler::OPERATION_UNARY_TRANS_TYPE;

//            p->rhs_.array_ = array;
//            p->rhs_.index_ = get_new_key(node.rhs_type_family_, key.first, node.rhs_.node_index_, RHS_NODE_TYPE);
//            p->rhs_.mapping_ = &mapping_;
//          }

//          void matrix_product_prototype(index_info const & key, statement_node const & node,  statement::container_type const * array) const{
//            mapped_matrix_product * p = new mapped_matrix_product("float");
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            p->lhs_.array_ = array;
//            p->lhs_.index_ = get_new_key(node.lhs_type_family_, key.first, node.lhs_.node_index_, LHS_NODE_TYPE);
//            p->lhs_.mapping_ = &mapping_;
////            p->is_lhs_transposed_ = array->at(node.lhs_.node_index_).op_type_ == scheduler::OPERATION_UNARY_TRANS_TYPE;

//            p->rhs_.array_ = array;
//            p->rhs_.index_ = get_new_key(node.rhs_type_family_, key.first, node.rhs_.node_index_, RHS_NODE_TYPE);
//            p->rhs_.mapping_ = &mapping_;
//          }

//          void scalar_reduction_prototype(index_info const & key, statement_node const & node,  statement::container_type const * array) const{
//            mapped_scalar_reduction * p = new mapped_scalar_reduction("float");
//            mapping_.insert(std::make_pair(key, tools::shared_ptr<mapped_container>(p)));
//            p->lhs_.array_ = array;
//            p->lhs_.index_ = get_new_key(node.lhs_type_family_, key.first, node.lhs_.node_index_, LHS_NODE_TYPE);
//            p->lhs_.mapping_ = &mapping_;
////            p->is_lhs_transposed_ = array->at(node.lhs_.node_index_).op_type_ == scheduler::OPERATION_UNARY_TRANS_TYPE;

//            p->rhs_.array_ = array;
//            p->rhs_.index_ = get_new_key(node.rhs_type_family_, key.first, node.rhs_.node_index_, RHS_NODE_TYPE);
//            p->rhs_.mapping_ = &mapping_;
//          }

//        public:

//          map_generate_prototype(std::map<void*, std::size_t> & memory, mapping_type & mapping, std::string & str, std::size_t & current_arg) : memory_(memory), mapping_(mapping), str_(str), current_arg_(current_arg) { }

//          void call_on_leaf(index_info const & key, statement_node const & node,  statement::container_type const * array) const {
//            if(key.second==PARENT_TYPE){
//              switch(node.op_type_){
//                case OPERATION_BINARY_ACCESS:
//                {
//                  index_info new_key = std::make_pair(key.first,LHS_NODE_TYPE);
//                  call_on_leaf(new_key, node, array);
//                  mapped_vector * v = static_cast<mapped_vector *>(mapping_.at(new_key).get());
//                  v->access_node_.array_ = array;
//                  v->access_node_.index_ = get_new_key(node.rhs_type_family_, key.first, node.rhs_.node_index_, RHS_NODE_TYPE);
//                  v->access_node_.mapping_ = &mapping_;
//                  break;
//                }
//                case OPERATION_BINARY_INNER_PROD_TYPE:
//                  scalar_reduction_prototype(key, node, array);
//                  break;
//                case OPERATION_BINARY_MAT_VEC_PROD_TYPE:
//                  vector_reduction_prototype(key, node, array);
//                  break;
//                case OPERATION_BINARY_MAT_MAT_PROD_TYPE:
//                  matrix_product_prototype(key, node, array);
//                  break;
//                default :
//                  throw "not handled";
//              }
//            }
//            else{
//              statement_node_type type;
//              statement_node_type_family type_family;
//              lhs_rhs_element element;
//              if(key.second==LHS_NODE_TYPE){
//                type = node.lhs_type_;
//                type_family = node.lhs_type_family_;
//                element = node.lhs_;
//              }
//              else{
//                type = node.rhs_type_;
//                type_family = node.rhs_type_family_;
//                element = node.rhs_;
//              }
//              switch(type_family){
//                case HOST_SCALAR_TYPE_FAMILY:
//                {
//                  switch(type){
//                    case HOST_SCALAR_FLOAT_TYPE : host_scalar_prototype(key, &element.host_float_); break;
//                    default : throw "not implemented";
//                  }
//                  break;
//                }
//                case SCALAR_TYPE_FAMILY:
//                {
//                  switch(type){
//                    case SCALAR_FLOAT_TYPE : scalar_prototype(key, element.scalar_float_); break;
//                    default : throw "not implemented";
//                  }
//                  break;
//                }
//                case VECTOR_TYPE_FAMILY:
//                {
//                  switch(type){
//                    case VECTOR_FLOAT_TYPE : vector_prototype(key, element.vector_float_);  break;
//                    default : throw "not implemented";
//                  }
//                  break;
//                }
//                case SYMBOLIC_VECTOR_TYPE_FAMILY:
//                {
//                  switch(type){
//                    case SYMBOLIC_VECTOR_FLOAT_TYPE : symbolic_vector_prototype(key, element.symbolic_vector_float_); break;
//                    default : throw "not implemented";
//                  }
//                  break;
//                }
//                case MATRIX_ROW_TYPE_FAMILY:
//                {
//                  switch(type){
//                    case MATRIX_ROW_FLOAT_TYPE : matrix_prototype(key, element.matrix_row_float_); break;
//                    default : throw "not implemented";
//                  }
//                  break;
//                }
//                case MATRIX_COL_TYPE_FAMILY:
//                {
//                  switch(type){
//                    case MATRIX_COL_FLOAT_TYPE : matrix_prototype(key, element.matrix_col_float_); break;
//                    default : throw "not implemented";
//                  }
//                  break;
//                }
//                case SYMBOLIC_MATRIX_TYPE_FAMILY:
//                {
//                  switch(type){
//                    case SYMBOLIC_MATRIX_FLOAT_TYPE : symbolic_matrix_prototype(key, element.symbolic_matrix_float_); break;
//                    default : throw "not implemented";
//                  }
//                  break;
//                }
//                default : throw "not implemented";
//              }
//           }
//        }

      };

    }

  }

}
#endif
