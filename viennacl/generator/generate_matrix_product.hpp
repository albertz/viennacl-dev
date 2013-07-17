#ifndef VIENNACL_GENERATOR_GENERATE_MATRIX_PRODUCT_HPP
#define VIENNACL_GENERATOR_GENERATE_MATRIX_PRODUCT_HPP

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


/** @file viennacl/generator/templates/matrix_product.hpp
 *
 * Kernel template for the vector reduction operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/detail.hpp"
#include "viennacl/generator/utils.hpp"

#include "viennacl/generator/template_base.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    class matrix_product : public template_base{
        typedef template_base base_type;
      public:
        typedef base_type::statements_type statements_type;

        class profile : public template_base::profile{
            friend class matrix_product;
            std::size_t lmem_used(std::size_t scalartype_size) const {
              std::size_t lmem_used = 0;
              if(use_LHS_shared_) lmem_used += (ml_ + 1) * (kl_ + 1) * scalartype_size;
              if(use_RHS_shared_) lmem_used += (kl_ + 1) * (nl_ + 1) * scalartype_size;
              return lmem_used;
            }

          public:
            /** @brief The user constructor */
            profile(unsigned int vectorization, unsigned int ml, unsigned int kl, unsigned int nl
                                               , unsigned int ms, unsigned int ks, unsigned int ns
                                               , bool use_LHS_shared, bool use_RHS_shared
                                               , unsigned int unroll) : template_base::profile(vectorization){
              ml_= ml; kl_=kl ; nl_=nl;
              ms_ = ms; ks_=ks; ns_=ns;
              use_LHS_shared_ = use_LHS_shared ; use_RHS_shared_ = use_RHS_shared;
              vectorization_ = vectorization;
              unroll_ = unroll;
            }

            void kernel_arguments(std::string & arguments_string) const{
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "M");
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "K");
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }

          private:
            unsigned int ml_;
            unsigned int kl_;
            unsigned int nl_;

            unsigned int ms_;
            unsigned int ks_;
            unsigned int ns_;

            bool use_LHS_shared_;
            bool use_RHS_shared_;

            unsigned int unroll_;
        };

      public:
        matrix_product( const & s, profile const & p) : template_base(s, profile_), profile_(p){ }

        void core(utils::kernel_generation_stream& stream) const{

        }

      private:
        profile profile_;
    };

}

#endif
