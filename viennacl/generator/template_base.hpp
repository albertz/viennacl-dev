#ifndef VIENNACL_GENERATOR_CODE_PROFILE_BASE
#define VIENNACL_GENERATOR_CODE_PROFILE_BASE

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


/** @file viennacl/generator/templates/profile_base.hpp
 *
 * Base classes for the profile
*/

#include <list>
#include <set>

#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/infos.hpp"

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/detail.hpp"

namespace viennacl{

  namespace generator{


    /** @brief Base class for an optimization profile */
    class template_base{
      public:
        typedef std::vector<scheduler::statement> statements_type;

        class profile{
          protected:
            virtual bool invalid_impl(viennacl::ocl::device const & dev, size_t scalartype_size) const{ return false; }
            virtual std::size_t lmem_used(std::size_t scalartype_size) const { return 0; }
          public:
            profile(unsigned int vectorization) : vectorization_(vectorization){ }

            virtual void kernel_arguments(std::string & arguments_string) const = 0;
            /** @brief returns whether or not the profile leads to undefined behavior on particular device
         *  @param dev the given device*/

            bool is_invalid(viennacl::ocl::device const & dev, size_t scalartype_size) const{
              //Query profile informations
              std::size_t size1, size2;
              set_local_sizes(size1, size2);

              //Query device informations
              size_t lmem_available = dev.local_mem_size();
              size_t max_workgroup_size = dev.max_work_group_size();
              std::vector<size_t> max_work_item_sizes = dev.max_work_item_sizes();

              bool invalid_work_group_sizes = size1*size2 > max_workgroup_size; // uses too much resources
              invalid_work_group_sizes = invalid_work_group_sizes || size1 > max_work_item_sizes[0];
              if(max_work_item_sizes.size()>1) invalid_work_group_sizes = invalid_work_group_sizes || size2 > max_work_item_sizes[1];

              return  invalid_work_group_sizes || lmem_used(scalartype_size)>lmem_available || invalid_impl(dev, scalartype_size);
            }
          protected:
            unsigned int vectorization_;
        };

      protected:
        statements_type const & statements_;
        mutable std::vector<detail::mapping_type> mapping_;
        profile const & profile_;

        template_base(statements_type const & s, profile const & p) : statements_(s), mapping_(s.size()), profile_(p) { }

      public:
        void prototype(utils::kernel_generation_stream & stream) const {
          std::map<void *, std::size_t> memory;
          std::string prototype;
          profile_.kernel_arguments(prototype);
          std::size_t current_arg = 0;
          std::size_t i = 0;
          for(statements_type::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it)
            detail::traverse(it->array(), detail::prototype_generation_traversal(memory, mapping_[i++], prototype, current_arg),false);
          prototype.erase(prototype.size()-1); //Last comma pruned
          stream << prototype << std::endl;
        }

    };


  }

}

#endif
