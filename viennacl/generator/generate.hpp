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
#include "viennacl/generator/generate_saxpy_vector.hpp"

namespace viennacl{

  namespace generator{


    using namespace viennacl::scheduler;

    enum operation_type_family{
      SAXPY_VECTOR,
      SAXPY_MATRIX,
      DOT,
      GEMV,
      GEMM
    };

    operation_type_family type_family_of(typename statement::container_type const & expr){
      switch (expr[0].lhs_type_family_) {
        case VECTOR_TYPE_FAMILY:
            return SAXPY_VECTOR;
          break;
        case MATRIX_ROW_TYPE_FAMILY:
          throw "not implemented";
          break;
        case MATRIX_COL_TYPE_FAMILY:
          throw "not implemented";
          break;
        default:
          throw "not implemented";
          break;
      }
    }

    class code_generator{
      private:
        typedef std::vector<scheduler::statement> statements_type;
      private:
        std::vector<detail::mapping_type> mapping_;
        statements_type statements_;
        saxpy_vector_profile saxpy_profile_;

      private:
        template<class Profile>
        void generate(Profile const & profile, utils::kernel_generation_stream & kss){
          //prototype:
          std::map<void *, std::size_t> memory;
          kss << "__kernel void kernel_0(" << std::endl;
          std::string prototype;
          profile.kernel_arguments(prototype);
          std::size_t current_arg = 0;
          std::size_t i = 0;
          for(typename statements_type::iterator it = statements_.begin() ; it != statements_.end() ; ++it)
            detail::traverse(it->array(), detail::prototype_generation_traversal(memory, mapping_[i++], prototype, current_arg),false);
          prototype.erase(prototype.size()-1); //Last comma pruned
          kss << prototype << std::endl;
          kss << ")" << std::endl;

          //core:
          kss << "{" << std::endl;
          kss.inc_tab();
          generate_saxpy_vector(profile, kss, statements_.begin(), statements_.end(), mapping_);
          kss.dec_tab();
          kss << "}" << std::endl;
        }

      public:
        code_generator() : saxpy_profile_(1,4,128) { }

        void add_statement(scheduler::statement const & s) { statements_.push_back(s); }

        std::string make_program_name(){
          std::string res;
          for(std::vector<scheduler::statement>::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it)
            detail::traverse(it->array(), detail::name_generation_traversal(res));
          return res;
        }

        std::string make_program_string(){
          std::size_t size = statements_.size();
          mapping_.resize(size);

          utils::kernel_generation_stream kss;

          //Headers generation
          kss << "#if defined(cl_khr_fp64)\n";
          kss <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
          kss <<  "#elif defined(cl_amd_fp64)\n";
          kss <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
          kss <<  "#endif\n";

          kss << std::endl;

          statement first_statement = statements_.front();
          switch(type_family_of(first_statement.array())){
            case SAXPY_VECTOR:  generate(saxpy_profile_, kss); break;
            default:  throw "not implemented";  break;
          }
          return kss.str();
        }
    };




  }

}
#endif
