#ifndef VIENNACL_GENERATOR_GENERATE_HPP
#define VIENNACL_GENERATOR_GENERATE_HPP

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


/** @file viennacl/generator/generate.hpp
    @brief User interface
*/

#include <cstring>
#include <vector>
#include <typeinfo>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/enqueue_tree.hpp"

#include "viennacl/generator/map_generate_prototype.hpp"

#include "viennacl/generator/generate_saxpy.hpp"
#include "viennacl/generator/generate_scalar_reduction.hpp"
#include "viennacl/generator/generate_vector_reduction.hpp"
#include "viennacl/generator/generate_matrix_product.hpp"

namespace viennacl{

  namespace generator{

    using namespace scheduler;

    enum operation_type_family{
      SCALAR_SAXPY,
      VECTOR_SAXPY,
      MATRIX_SAXPY,
      SCALAR_REDUCE,
      VECTOR_REDUCE,
      MATRIX_PRODUCT
    };

    unsigned int count(typename statement::container_type const & expr, operation_node_type op_type){
      unsigned int res = 0;
      for(typename statement::container_type::const_iterator it = expr.begin() ; it != expr.end() ; ++it)
        res += static_cast<unsigned int>(it->op_type_==op_type);
      return res;
    }

    operation_type_family type_family_of(typename statement::container_type const & expr){
      switch(expr[0].lhs_type_family_){
        case VECTOR_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_MAT_VEC_PROD_TYPE))
            return VECTOR_REDUCE;
          else
            return VECTOR_SAXPY;
        case MATRIX_ROW_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_MAT_MAT_PROD_TYPE))
            return MATRIX_PRODUCT;
          else
            return MATRIX_SAXPY;
        case MATRIX_COL_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_MAT_MAT_PROD_TYPE))
            return MATRIX_PRODUCT;
          else
            return MATRIX_SAXPY;
        case SCALAR_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_INNER_PROD_TYPE))
            return SCALAR_REDUCE;
          else
            return SCALAR_SAXPY;
        default:
          throw "not implemented";
      }
    }

    class code_generator{
      private:
        typedef std::list<std::pair<operation_type_family, generator::template_base::statements_type> > statements_type;

        template<class T>
        static void merge(T & first, T const & second){
          first.insert(first.end(), second.begin(), second.end());
        }

        struct vec_repr{
            typedef const char * result_type;
            template<class ScalarType>
            result_type operator()(vector_base<ScalarType> const & t) const {
              if(t.start()>0)
                return "vr";
              if(t.stride()>0)
                return "vst";
              return "v";
            }
        };

        struct mat_repr{
            typedef const char * result_type;
            template<class ScalarType, class Layout>
            result_type operator()(matrix_base<ScalarType, Layout> const & t) const {
              if(t.start1()>0)
                return "mr";
              return "m";
            }
        };

      public:
        code_generator() : vector_saxpy_profile_(1,128,128,true)
                          , matrix_saxpy_profile_(1,16,16,16,16,true)
                          , scalar_reduction_profile_(1, 128, 128, true)
                          , vector_reduction_profile_(1, 1, 256, 32)
                          , matrix_product_profile_(1,32,32,32,4,4,4,true,false,1)
                           { }

        bool add_statement(scheduler::statement const & s) {
          operation_type_family expr_type = type_family_of(s.array());
          if(statements_.empty())
            statements_.push_back(std::make_pair(expr_type,template_base::statements_type(1,s)));
          else
            if(statements_.back().first == expr_type)
              statements_.back().second.push_back(s);
            else
              statements_.push_back(std::make_pair(expr_type,template_base::statements_type(1,s)));
          return true;
        }

        void configure_program(viennacl::ocl::program & p){
          unsigned int kernel_id = 0;
          const char * kernel_prefix = "kernel_";
          for(statements_type::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it){
            switch(it->first){
              case VECTOR_SAXPY:
                vector_saxpy_profile_.enqueue_kernel_arguments(it->second,p.get_kernel(kernel_prefix+utils::to_string(kernel_id++)));
                break;
              case MATRIX_SAXPY:
                break;
              case SCALAR_REDUCE:
                break;
              case VECTOR_REDUCE:
                break;
              case MATRIX_PRODUCT:
                break;
              default:
                break;
            }
          }
        }

        void make_program_name(char *& ptr) const{
          for(statements_type::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it){
            for(std::vector<scheduler::statement>::const_iterator iit = it->second.begin() ; iit != it->second.end() ; ++iit){
              scheduler::statement::container_type const & expr = iit->array();
              for(std::size_t j = 0 ; j < expr.size() ; ++j){
                if(expr[j].lhs_type_family_==scheduler::COMPOSITE_OPERATION_FAMILY){
                  strcat(ptr, "_");
                }
                else{
                  if(expr[j].lhs_type_family_==scheduler::VECTOR_TYPE_FAMILY)
                    strcat(ptr, utils::call_on_vector(expr[j].lhs_type_, expr[j].lhs_, vec_repr()));
                  else if(expr[j].lhs_type_family_==scheduler::MATRIX_ROW_TYPE_FAMILY){
                    strcat(ptr, utils::call_on_matrix(expr[j].lhs_type_, expr[j].lhs_, mat_repr()));
                  }
                  strcat(ptr, detail::generate_scalartype(expr[j].lhs_type_));
                }

                strcat(ptr, detail::generate(expr[j].op_type_));

                if(expr[j].rhs_type_family_==scheduler::COMPOSITE_OPERATION_FAMILY){
                  strcat(ptr, "_");
                }
                else{
                  if(expr[j].rhs_type_family_==scheduler::VECTOR_TYPE_FAMILY)
                    strcat(ptr, utils::call_on_vector(expr[j].rhs_type_, expr[j].rhs_, vec_repr()));
                  else if(expr[j].rhs_type_family_==scheduler::MATRIX_ROW_TYPE_FAMILY){
                    strcat(ptr, utils::call_on_matrix(expr[j].rhs_type_, expr[j].rhs_, mat_repr()));
                  }
                  strcat(ptr, detail::generate_scalartype(expr[j].rhs_type_));
                }
              }
            }
          }
        }



        std::string make_program_string(std::vector<std::string> & kernels_name){
          unsigned int kernel_name_offset = 0;
          utils::kernel_generation_stream stream;

          //Headers generation
          stream << "#if defined(cl_khr_fp64)\n";
          stream <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
          stream <<  "#elif defined(cl_amd_fp64)\n";
          stream <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
          stream <<  "#endif\n";
          stream << std::endl;

          for(statements_type::iterator it = statements_.begin() ; it != statements_.end() ; ++it){
            switch(it->first){
              case VECTOR_SAXPY:
                merge(kernels_name, vector_saxpy(it->second, vector_saxpy_profile_)(stream, kernel_name_offset));
                break;
              case MATRIX_SAXPY:
                merge(kernels_name, matrix_saxpy(it->second, matrix_saxpy_profile_)(stream, kernel_name_offset));
                break;
              case SCALAR_REDUCE:
                merge(kernels_name, scalar_reduction(it->second, scalar_reduction_profile_)(stream, kernel_name_offset));
                break;
              case VECTOR_REDUCE:
                merge(kernels_name, vector_reduction(it->second, vector_reduction_profile_)(stream, kernel_name_offset));
                break;
              case MATRIX_PRODUCT:
                merge(kernels_name, matrix_product(it->second,matrix_product_profile_)(stream, kernel_name_offset));
                break;
              default:
                break;
            }
          }
          return stream.str();
        }

      private:
        statements_type statements_;

        vector_saxpy::profile vector_saxpy_profile_;
        matrix_saxpy::profile matrix_saxpy_profile_;
        scalar_reduction::profile scalar_reduction_profile_;
        vector_reduction::profile vector_reduction_profile_;
        matrix_product::profile matrix_product_profile_;

    };




  }

}
#endif
