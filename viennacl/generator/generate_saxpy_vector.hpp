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

    void generate_saxpy_vector(saxpy_vector_profile const & prof, utils::kernel_generation_stream& kss, std::vector<viennacl::scheduler::statement> const & statements, std::vector< std::vector<detail::symbolic_container> > const & symbolic_mappings){

      kss << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0))" << std::endl;
      kss << "{" << std::endl;
      kss.inc_tab();

      for(unsigned int i = 0 ; i < statements.size() ; ++i){
        detail::traverse(statements[i].array(), 0, detail::fetching_traversal(kss,symbolic_mappings[i],"i"));
      }
//      //Set access indices
//      for(std::list<tools::shared_ptr<symbolic_binary_expression_tree_infos_base> >::iterator it=expressions_.begin() ; it!=expressions_.end();++it){
//        (*it)->access_index("i","0");
//        (*it)->fetch(kss);
//      }

//      //Compute expressions
//      for(std::list<tools::shared_ptr<symbolic_binary_expression_tree_infos_base> >::iterator it=expressions_.begin() ; it!=expressions_.end();++it)
//        kss << (*it)->generate() << ";" << std::endl;
//      for(std::list<tools::shared_ptr<symbolic_binary_expression_tree_infos_base> >::iterator it=expressions_.begin() ; it!=expressions_.end();++it)
//        (*it)->write_back(kss);

      kss.dec_tab();
      kss << "}" << std::endl;

//      for(std::list<tools::shared_ptr<symbolic_binary_expression_tree_infos_base> >::iterator it = expressions_.begin(); it != expressions_.end() ; ++it)
//        (*it)->clear_private_value();
    }

  }

}

#endif
