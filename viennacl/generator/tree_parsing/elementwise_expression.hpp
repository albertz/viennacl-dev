#ifndef VIENNACL_GENERATOR_GENERATE_EXPRESSION_GENERATION_HPP
#define VIENNACL_GENERATOR_GENERATE_EXPRESSION_GENERATION_HPP

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


/** @file viennacl/generator/helpers.hpp
    @brief several code generation helpers
*/

#include <set>

#include "CL/cl.h"

#include "viennacl/forwards.h"
#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/forwards.h"
#include "viennacl/generator/tree_parsing/traverse.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{

      /** @brief generate a string from an operation_node_type */
      static const char * generate(operation_node_type type){
        // unary expression
        switch(type){
          //Function
          case OPERATION_UNARY_ABS_TYPE : return "abs";
          case OPERATION_BINARY_ELEMENT_POW_TYPE : return "pow";

          //Arithmetic
          case OPERATION_BINARY_ASSIGN_TYPE : return "=";
          case OPERATION_BINARY_INPLACE_ADD_TYPE : return "+=";
          case OPERATION_BINARY_INPLACE_SUB_TYPE : return "-=";
          case OPERATION_BINARY_ADD_TYPE : return "+";
          case OPERATION_BINARY_SUB_TYPE : return "-";
          case OPERATION_BINARY_MULT_TYPE : return "*";
          case OPERATION_BINARY_DIV_TYPE : return "/";
          case OPERATION_BINARY_ACCESS_TYPE : return "[]";

          //Relational
          case OPERATION_BINARY_ELEMENT_EQ_TYPE : return "isequal";
          case OPERATION_BINARY_ELEMENT_NEQ_TYPE : return "isnotequal";
          case OPERATION_BINARY_ELEMENT_GREATER_TYPE : return "isgreater";
          case OPERATION_BINARY_ELEMENT_GEQ_TYPE : return "isgreaterequal";
          case OPERATION_BINARY_ELEMENT_LESS_TYPE : return "isless";
          case OPERATION_BINARY_ELEMENT_LEQ_TYPE : return "islessequal";

          case OPERATION_BINARY_ELEMENT_FMAX_TYPE : return "fmax";
          case OPERATION_BINARY_ELEMENT_FMIN_TYPE : return "fmin";

          //Unary
          case OPERATION_BINARY_REDUCE_TYPE : return "reduce";

          case OPERATION_UNARY_TRANS_TYPE : return "trans";
          case OPERATION_UNARY_ROW_WISE_TYPE : return "rowwise";
          case OPERATION_UNARY_COLUMN_WISE_TYPE : return "columnwise";

          //Binary
          //Leaves
          case OPERATION_BINARY_INNER_PROD_TYPE : return "iprod";
          case OPERATION_BINARY_MAT_MAT_PROD_TYPE : return "mmprod";
          case OPERATION_BINARY_MAT_VEC_PROD_TYPE : return "mvprod";


          default : throw generator_not_supported_exception("Unsupported operator");
        }
      }


      /** @brief functor for generating the expression string from a statement */
      class expression_generation_traversal : public traversal_functor{
        private:
          std::pair<std::string, std::string> index_string_;
          int simd_element_;
          std::string & str_;
          mapping_type const & mapping_;

        public:
          expression_generation_traversal(std::pair<std::string, std::string> const & index, int simd_element, std::string & str, mapping_type const & mapping) : index_string_(index), simd_element_(simd_element), str_(str), mapping_(mapping){ }

          void call_before_expansion(scheduler::statement_node const * root_node) const
          {
              if(root_node->op.type_subfamily==scheduler::OPERATION_ELEMENTWISE_FUNCTION_TYPE_SUBFAMILY)
                  str_+=generate(root_node->op.type);
              str_+="(";
          }
          void call_after_expansion(scheduler::statement_node const *) const
          { str_+=")";
          }

          void operator()(scheduler::statement const *, scheduler::statement_node const * root_node, detail::node_type node_type) const {
            if(node_type==PARENT_NODE_TYPE)
            {
               switch(root_node->op.type_subfamily){
                   case scheduler::OPERATION_STRUCTUREWISE_FUNCTION_TYPE_SUBFAMILY: str_ += generate(index_string_, simd_element_, *mapping_.at(std::make_pair(root_node, node_type))); break;
                   case scheduler::OPERATION_ELEMENTWISE_OPERATOR_TYPE_SUBFAMILY: str_ += generate(root_node->op.type); break;
                   case scheduler::OPERATION_ELEMENTWISE_FUNCTION_TYPE_SUBFAMILY: str_ += ","; break;
                   default: break;
               }
            }
            else{
              if(node_type==LHS_NODE_TYPE){
                if(root_node->lhs.type_family!=scheduler::COMPOSITE_OPERATION_FAMILY)
                  str_ += detail::generate(index_string_,simd_element_, *mapping_.at(std::make_pair(root_node,node_type)));
              }
              else if(node_type==RHS_NODE_TYPE){
                if(root_node->rhs.type_family!=scheduler::COMPOSITE_OPERATION_FAMILY)
                  str_ += detail::generate(index_string_,simd_element_, *mapping_.at(std::make_pair(root_node,node_type)));
              }
            }
          }
      };

      static void generate_all_lhs(scheduler::statement const & statement
                                , scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , int simd_element
                                , std::string & str
                                , detail::mapping_type const & mapping){
        if(root_node.lhs.type_family==scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.lhs.node_index], detail::expression_generation_traversal(index, simd_element, str, mapping));
        else
          str += detail::generate(index, simd_element,*mapping.at(std::make_pair(&root_node,detail::LHS_NODE_TYPE)));
      }


      static void generate_all_rhs(scheduler::statement const & statement
                                , scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , int simd_element
                                , std::string & str
                                , detail::mapping_type const & mapping){
        if(root_node.rhs.type_family==scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.rhs.node_index], detail::expression_generation_traversal(index, simd_element, str, mapping));
        else
          str += detail::generate(index, simd_element,*mapping.at(std::make_pair(&root_node,detail::RHS_NODE_TYPE)));
      }



    }
  }
}
#endif
