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

#include <vector>
#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/generate_saxpy.hpp"
#include "viennacl/generator/generate_scalar_reduction.hpp"
#include "viennacl/generator/generate_vector_reduction.hpp"

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
    }

    operation_type_family type_family_of(typename statement::container_type const & expr){
      switch(expr[0].lhs_type_family_){
        case VECTOR_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_PROD_TYPE))
            return VECTOR_REDUCE;
          else
            return VECTOR_SAXPY;
        case MATRIX_ROW_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_PROD_TYPE))
            return MATRIX_PRODUCT;
          else
            return MATRIX_SAXPY;
        case MATRIX_COL_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_PROD_TYPE))
            return MATRIX_PRODUCT;
          else
            return MATRIX_SAXPY;
        case SCALAR_TYPE_FAMILY :
          if(count(expr, OPERATION_BINARY_INNER_PROD_TYPE))
            return SCALAR_REDUCE;
          else
            return SCALAR_SAXPY;
        default: throw "not implemented";
      }
    }

    class code_generator{
      private:
        generator::template_base::statements_type statements_;
        vector_saxpy::profile vector_saxpy_profile_;
        matrix_saxpy::profile matrix_saxpy_profile_;
        vector_reduction::profile vector_reduction_profile_;


      private:
        template<class Generator>
        static void generate(Generator const & g, utils::kernel_generation_stream & stream){

          //prototype:
          stream << "__kernel void kernel_0(" << std::endl;
          g.prototype(stream);
          stream << ")" << std::endl;

          //core:
          stream << "{" << std::endl;
          stream.inc_tab();
          g.core(stream);
          stream.dec_tab();
          stream << "}" << std::endl;
        }

      public:
        code_generator() : vector_saxpy_profile_(1,128,128,true), matrix_saxpy_profile_(1,16,16,16,16,true), vector_reduction_profile_(1, 1, 256, 32) { }

        void add_statement(scheduler::statement const & s) { statements_.push_back(s); }

        std::string make_program_name(){
          return "";
        }

        std::string make_program_string(){
          utils::kernel_generation_stream stream;

          //Headers generation
          stream << "#if defined(cl_khr_fp64)\n";
          stream <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
          stream <<  "#elif defined(cl_amd_fp64)\n";
          stream <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
          stream <<  "#endif\n";

          stream << std::endl;

          statement first_statement = statements_.front();
          switch(type_family_of(first_statement.array())){
            case VECTOR_SAXPY:  generate(vector_saxpy(statements_, vector_saxpy_profile_), stream); break;
            case MATRIX_SAXPY:  generate(matrix_saxpy(statements_, matrix_saxpy_profile_), stream); break;
            case VECTOR_REDUCE: generate(vector_reduction(statements_, vector_reduction_profile_), stream); break;
            default:  throw "not implemented";  break;
          }
          return stream.str();
        }
    };




  }

}
#endif
