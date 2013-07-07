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
namespace viennacl{

  namespace generator{

    using namespace viennacl::scheduler;

    std::string make_program_name(scheduler::statement const & statement){
      std::string res;
      detail::traverse(statement.array(), 0, detail::generation_traversal(res));
      return res;
    }

    std::string make_program_string(scheduler::statement const & statement){
      std::map<cl_mem, unsigned int> memory;
      std::string str;
      detail::traverse(statement.array(), 0, detail::header_generation_traversal(memory, str));
      str.erase(str.size()-1);
      return str;
    }

  }

}
#endif
