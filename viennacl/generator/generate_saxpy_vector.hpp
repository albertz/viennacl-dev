#ifndef VIENNACL_GENERATOR_GENERATE_SAXPY_VECTOR_HPP
#define VIENNACL_GENERATOR_GENERATE_SAXPY_VECTOR_HPP

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


/** @file viennacl/generator/templates/saxpy.hpp
 *
 * Kernel template for the SAXPY operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/profile_base.hpp"
#include "viennacl/generator/detail.hpp"
#include "viennacl/generator/utils.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    class saxpy_vector_profile : public profile_base{
      public:
        /** @brief The constructor */
        saxpy_vector_profile(unsigned int vectorization, unsigned int loop_unroll, size_t group_size) : profile_base(vectorization){
          loop_unroll_ = loop_unroll;
          group_size_ = group_size;
        }

        void set_local_sizes(std::size_t & x, std::size_t & y) const{
          x = group_size_;
          y = 1;
        }

        static void kernel_arguments(std::string & arguments_string){
          arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
        }
      private:
        unsigned int loop_unroll_;
        unsigned int group_size_;
    };

    template<class InputIterator>
    void generate_saxpy_vector(saxpy_vector_profile const & prof, utils::kernel_generation_stream& stream, InputIterator first, InputIterator last, std::vector<detail::mapping_type> & mapping){

      stream << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0))" << std::endl;
      stream << "{" << std::endl;
      stream.inc_tab();

      //Fetches entries to registers
      std::set<std::string>  fetched;
      for(std::vector<detail::mapping_type>::iterator it = mapping.begin() ; it != mapping.end() ; ++it)
        for(detail::mapping_type::reverse_iterator it2 = it->rbegin() ; it2 != it->rend() ; ++it2)
          it2->second->fetch(fetched, stream);


      std::size_t i = 0;
      for(InputIterator it = first ; it != last ; ++it){
          detail::traverse(it->array(), detail::expression_generation_traversal(stream,mapping[i++]), true);
          stream << ";" << std::endl;
      }

      //Writes back
      for(std::vector<detail::mapping_type>::iterator it = mapping.begin() ; it != mapping.end() ; ++it)
        for(detail::mapping_type::iterator it2 = it->begin() ; it2 != it->end() ; ++it2)
          it2->second->write_back(fetched, stream);

      stream.dec_tab();
      stream << "}" << std::endl;
    }

  }

}

#endif
