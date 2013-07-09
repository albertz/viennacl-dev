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

    std::string make_program_name(std::vector<scheduler::statement> const & statements){
      std::string res;
      for(std::vector<scheduler::statement>::const_iterator it = statements.begin() ; it != statements.end() ; ++it)
        detail::traverse(it->array(), 0, detail::name_generation_traversal(res));
      return res;
    }

    std::string make_program_string(std::vector<scheduler::statement> const & statements){
      std::map<cl_mem, std::size_t> memory;
      std::size_t size = statements.size();
      std::vector<detail::mapping_type> mapping(size);

      utils::kernel_generation_stream kss;

      //Headers generation
      kss << "#if defined(cl_khr_fp64)\n";
      kss <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
      kss <<  "#elif defined(cl_amd_fp64)\n";
      kss <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
      kss <<  "#endif\n";

      kss << std::endl;

      //Prototype generation
      kss << "__kernel void kernel_0(" << std::endl;
      std::string prototype;
      for(std::size_t i = 0 ; i < size ; ++i)
        detail::traverse(statements[i].array(), 0, detail::prototype_generation_traversal(memory, mapping[i], prototype));
      prototype.erase(prototype.size()-1); //Last comma pruned
      kss << prototype << std::endl;
      kss << ")" << std::endl;
      kss << "{" << std::endl;
      kss.inc_tab();

      //body generation
      {
        generate_saxpy_vector(saxpy_vector_profile(1,4,128), kss, statements, mapping);
      }

      kss.dec_tab();
      kss << "}" << std::endl;
      return kss.str();
    }

  }

}
#endif
