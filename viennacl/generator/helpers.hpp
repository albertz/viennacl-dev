#ifndef VIENNACL_GENERATOR_GENERATE_UTILS_HPP
#define VIENNACL_GENERATOR_GENERATE_UTILS_HPP

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

#include "viennacl/generator/utils.hpp"
#include "viennacl/generator/forwards.h"

namespace viennacl{

  namespace generator{

    namespace detail{

    /** @brief generate the string for a pointer kernel argument */
      static std::string generate_value_kernel_argument(std::string const & scalartype, std::string const & name){
        return scalartype + ' ' + name + ",";
      }

      /** @brief generate the string for a pointer kernel argument */
      static std::string generate_pointer_kernel_argument(std::string const & address_space, std::string const & scalartype, std::string const & name){
        return address_space +  " " + scalartype + "* " + name + ",";
      }

      /** @brief generate a string from an operation_node_type */
      static const char * generate(operation_node_type type){
        // unary expression
        switch(type){
          //Function
          case OPERATION_UNARY_ABS_TYPE : return "abs";
          case OPERATION_BINARY_ELEMENT_POW_TYPE : return "pow";

          //Unary
          case OPERATION_UNARY_TRANS_TYPE : return "trans";

          //Binary
          //Leaves
          case OPERATION_BINARY_INNER_PROD_TYPE : return "iprod";
          case OPERATION_BINARY_MAT_MAT_PROD_TYPE : return "mmprod";
          case OPERATION_BINARY_MAT_VEC_PROD_TYPE : return "mvprod";

          //Arithmetic
          case OPERATION_BINARY_ASSIGN_TYPE : return "=";
          case OPERATION_BINARY_INPLACE_ADD_TYPE : return "+=";
          case OPERATION_BINARY_INPLACE_SUB_TYPE : return "-=";
          case OPERATION_BINARY_ADD_TYPE : return "+";
          case OPERATION_BINARY_SUB_TYPE : return "-";
          case OPERATION_BINARY_MULT_TYPE : return "*";
          case OPERATION_BINARY_DIV_TYPE : return "/";
          case OPERATION_BINARY_ACCESS_TYPE : return "[]";

          //Binary elementwise
          case OPERATION_BINARY_ELEMENT_EQ_TYPE : return "==";
          case OPERATION_BINARY_ELEMENT_NEQ_TYPE : return "!=";
          case OPERATION_BINARY_ELEMENT_GREATER_TYPE : return ">";
          case OPERATION_BINARY_ELEMENT_GEQ_TYPE : return ">=";
          case OPERATION_BINARY_ELEMENT_LESS_TYPE : return "<";
          case OPERATION_BINARY_ELEMENT_LEQ_TYPE : return "<=";


          default : throw "not implemented";
        }
      }

      /** @brief Recursively execute a functor on a statement */
      template<class Fun>
      static void traverse(scheduler::statement const & statement, scheduler::statement_node const & root_node, Fun const & fun, bool recurse_structurewise_function /* see forwards.h for default argument */){

        if(root_node.op.type_family==OPERATION_UNARY_TYPE_FAMILY)
        {
          //Self:
          fun(&statement, &root_node, PARENT_NODE_TYPE);

          //Lhs:
          fun.call_before_expansion(&root_node);

          if(root_node.lhs.type_family==COMPOSITE_OPERATION_FAMILY)
              traverse(statement, statement.array()[root_node.lhs.node_index], fun, recurse_structurewise_function);
          fun(&statement, &root_node, LHS_NODE_TYPE);

          fun.call_after_expansion(&root_node);
        }
        else if(root_node.op.type_family==OPERATION_BINARY_TYPE_FAMILY)
        {
          bool deep_recursion = recurse_structurewise_function || root_node.op.type_subfamily!=scheduler::OPERATION_STRUCTUREWISE_FUNCTION_TYPE_SUBFAMILY;

          fun.call_before_expansion(&root_node);

          //Lhs:
          if(deep_recursion){
            if(root_node.lhs.type_family==COMPOSITE_OPERATION_FAMILY)
              traverse(statement, statement.array()[root_node.lhs.node_index], fun, recurse_structurewise_function);
            fun(&statement, &root_node, LHS_NODE_TYPE);
          }

          //Self:
          fun(&statement, &root_node, PARENT_NODE_TYPE);

          //Rhs:
          if(deep_recursion){
            if(root_node.rhs.type_family==COMPOSITE_OPERATION_FAMILY)
              traverse(statement, statement.array()[root_node.rhs.node_index], fun, recurse_structurewise_function);
            fun(&statement, &root_node, RHS_NODE_TYPE);
          }

          fun.call_after_expansion(&root_node);

        }
      }

      /** @brief base functor class for traversing a statement */
      class traversal_functor{
        public:
          void call_before_expansion(scheduler::statement_node const *) const { }
          void call_after_expansion(scheduler::statement_node const *) const { }
      };

      /** @brief functor for generating the prototype of a statement */
      class prototype_generation_traversal : public traversal_functor{
        private:
          std::string & str_;
          std::set<std::string> & already_generated_;
          mapping_type const & mapping_;
        public:
          prototype_generation_traversal(std::set<std::string> & already_generated, std::string & str, mapping_type const & mapping) : already_generated_(already_generated), str_(str),  mapping_(mapping){ }

          void operator()(scheduler::statement const *, scheduler::statement_node const * root_node, detail::node_type node_type) const {
              if( (node_type==detail::LHS_NODE_TYPE && root_node->lhs.type_family!=scheduler::COMPOSITE_OPERATION_FAMILY)
                ||(node_type==detail::RHS_NODE_TYPE && root_node->rhs.type_family!=scheduler::COMPOSITE_OPERATION_FAMILY) )
                  append_kernel_arguments(already_generated_, str_, *mapping_.at(std::make_pair(root_node,node_type)));
          }
      };

      /** @brief functor for fetching the elements of a statement */
      class fetch_traversal : public traversal_functor{
        private:
          std::set<std::string> & fetched_;
          std::pair<std::string, std::string> index_string_;
          utils::kernel_generation_stream & stream_;
          mapping_type const & mapping_;
        public:
          fetch_traversal(std::set<std::string> & fetched, std::pair<std::string, std::string> const & index, utils::kernel_generation_stream & stream, mapping_type const & mapping) : fetched_(fetched), index_string_(index), stream_(stream), mapping_(mapping){ }

          void operator()(scheduler::statement const *, scheduler::statement_node const * root_node, detail::node_type node_type) const {
            if( (node_type==detail::LHS_NODE_TYPE && root_node->lhs.type_family!=scheduler::COMPOSITE_OPERATION_FAMILY)
              ||(node_type==detail::RHS_NODE_TYPE && root_node->rhs.type_family!=scheduler::COMPOSITE_OPERATION_FAMILY) )
              fetch(index_string_, fetched_, stream_, *mapping_.at(std::make_pair(root_node, node_type)));
          }
      };

      /** @brief functor for fetching the LHS of a statement's node
      *
      *   Forwards to fetch_traversal functor if the LHS is not a leaf
      */
      static void fetch_all_lhs(std::set<std::string> & fetched
                                , scheduler::statement const & statement
                                , scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , utils::kernel_generation_stream & stream
                                , detail::mapping_type const & mapping){
        if(root_node.lhs.type_family==scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.lhs.node_index], detail::fetch_traversal(fetched, index, stream, mapping));
        else
          detail::fetch(index, fetched, stream, *mapping.at(std::make_pair(&root_node,detail::LHS_NODE_TYPE)));

      }

      /** @brief functor for fetching the RHS of a statement's node
      *
      *   Forwards to fetch_traversal functor if the RHS is not a leaf
      */
      static void fetch_all_rhs(std::set<std::string> & fetched
                                , scheduler::statement const & statement
                                , scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , utils::kernel_generation_stream & stream
                                , detail::mapping_type const & mapping){
        if(root_node.rhs.type_family==scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.rhs.node_index], detail::fetch_traversal(fetched, index, stream, mapping));
        else
          detail::fetch(index, fetched, stream, *mapping.at(std::make_pair(&root_node,detail::RHS_NODE_TYPE)));

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
