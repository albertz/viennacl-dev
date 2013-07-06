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

    class generate_functor{
        std::string & str_;
      public:
        generate_functor(std::string & str) : str_(str){ }
        template<class T> void operator()(T const & t) const{ str_+=detail::generate(t); }
    };

    std::string make_program_name(viennacl::scheduler::statement const & statement){
      std::string res;
      detail::for_each(statement,generate_functor(res));
      return res;
    }

    template<class ExpressionsContainer>
    std::string make_program_string(ExpressionsContainer const & statements){
      return "";
    }

  }

}
#endif
