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

#include "viennacl/generator/detail.hpp"
#include "viennacl/generator/saxpy.hpp"

namespace viennacl{

  namespace generator{


    using namespace viennacl::scheduler;

    enum operation_type_family{
      VECTOR_SAXPY,
      MATRIX_SAXPY,
      DOT,
      GEMV,
      GEMM
    };

    operation_type_family type_family_of(typename statement::container_type const & expr){
      switch(expr[0].lhs_type_family_){
        case VECTOR_TYPE_FAMILY : return VECTOR_SAXPY;
        case MATRIX_ROW_TYPE_FAMILY : return MATRIX_SAXPY;
        case MATRIX_COL_TYPE_FAMILY : return MATRIX_SAXPY;
        case SCALAR_TYPE_FAMILY : return DOT;
        default: throw "not implemented";
      }
    }

    class code_generator{
      private:
        template_base::statements_type statements_;
        vector_saxpy::profile vector_saxpy_profile_;
        matrix_saxpy::profile matrix_saxpy_profile_;


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
        code_generator() : vector_saxpy_profile_(1,128,128,true), matrix_saxpy_profile_(1,16,16,16,16,true) { }

        void add_statement(scheduler::statement const & s) { statements_.push_back(s); }

        std::string make_program_name(){
          std::string res;
          for(std::vector<scheduler::statement>::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it)
            detail::traverse(it->array(), detail::name_generation_traversal(res));
          return res;
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
            default:  throw "not implemented";  break;
          }
          return stream.str();
        }
    };




  }

}
#endif
