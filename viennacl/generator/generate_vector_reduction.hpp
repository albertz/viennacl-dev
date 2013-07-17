#ifndef VIENNACL_GENERATOR_GENERATE_VECTOR_REDUCTION_HPP
#define VIENNACL_GENERATOR_GENERATE_VECTOR_REDUCTION_HPP

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


/** @file viennacl/generator/templates/vector_reduction.hpp
 *
 * Kernel template for the vector reduction operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/detail.hpp"
#include "viennacl/generator/utils.hpp"

#include "viennacl/generator/generate_template_base.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    class vector_reduction : public template_base{
        typedef template_base base_type;
      public:
        typedef base_type::statements_type statements_type;

        class profile : public template_base::profile{
            friend class vector_reduction;
            std::size_t lmem_used(std::size_t scalartype_size) const {
              return m_*(k_+1)*scalartype_size;
            }

          public:
            /** @brief The user constructor */
            profile(unsigned int vectorization, unsigned int m, unsigned int k, unsigned int num_groups) : template_base::profile(vectorization), m_(m), k_(k), num_groups_(num_groups){ }
            void set_local_sizes(std::size_t& s1, std::size_t& s2) const{
              s1 = m_;
              s2 = k_;
            }
            void kernel_arguments(std::string & arguments_string) const{
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "M");
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }
          private:
            unsigned int m_;
            unsigned int k_;
            unsigned int num_groups_;
        };

      public:
        vector_reduction(template_base::statements_type const & s, profile const & p) : template_base(s, profile_), profile_(p){ }

        void core(utils::kernel_generation_stream& stream) const{

        }

      private:
        profile profile_;
    };
  }
}

#endif
