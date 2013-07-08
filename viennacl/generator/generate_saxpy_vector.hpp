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

        /** @brief The user constructor */
        saxpy_vector_profile(unsigned int vectorization, unsigned int loop_unroll, size_t group_size0) : profile_base(vectorization){
          loop_unroll_ = loop_unroll;
          group_size_ = group_size0;
        }

        /** @brief Returns the unrolling factor */
        unsigned int loop_unroll() const{ return loop_unroll_; }

        /** @brief Return the group sizes used by this kernel */
        std::pair<size_t,size_t> local_work_size() const{ return std::make_pair(group_size_,1); }

        /** @brief returns whether or not the profile leads to undefined behavior on particular device
         *  @param dev the given device*/
        bool is_invalid(viennacl::ocl::device const & dev, size_t scalartype_size) const {
          return profile_base::invalid_base(dev,0);
        }

      private:
        unsigned int loop_unroll_;
        unsigned int group_size_;
    };

    void generate_saxpy_vector(saxpy_vector_profile const & prof, utils::kernel_generation_stream& kss, std::vector<viennacl::scheduler::statement> const & statements, std::vector<detail::mapping_type> const & mapping){

      kss << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0))" << std::endl;
      kss << "{" << std::endl;
      kss.inc_tab();

      //Fetches entries to registers
      std::set<std::string> fetched;
      for(std::vector<detail::mapping_type>::const_iterator it = mapping.begin() ; it != mapping.end() ; ++it)
        for(detail::mapping_type::const_iterator it2 = it->begin() ; it2 != it->end() ; ++it2)
          if(fetched.insert(it2->second.name_).second)
            kss << it2->second.scalartype_ << " " << it2->second.name_ << "_private" << "=" << it2->second.name_ + "[i];" << std::endl;


      for(std::size_t i = 0 ; i < statements.size() ; ++i){
        detail::traverse(statements[i].array(),0,detail::expression_generation_traversal<detail::add_suffix_to_name>(kss,mapping[i],detail::add_suffix_to_name("_private")));
        kss << ";" << std::endl;
      }

      //Writes back
      for(std::set<std::string>::iterator it = fetched.begin() ; it != fetched.end() ; ++it)
        kss <<  *it + "[i]" << '=' << *it << "_private" << ';' << std::endl;

      kss.dec_tab();
      kss << "}" << std::endl;

//      for(std::list<tools::shared_ptr<symbolic_binary_expression_tree_infos_base> >::iterator it = expressions_.begin(); it != expressions_.end() ; ++it)
//        (*it)->clear_private_value();
    }

  }

}

#endif
