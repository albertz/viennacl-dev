#ifndef VIENNACL_GENERATOR_GENERATE_SCALAR_REDUCTION_HPP
#define VIENNACL_GENERATOR_GENERATE_SCALAR_REDUCTION_HPP

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


/** @file viennacl/generator/templates/scalar_reduction.hpp
 *
 * Kernel template for the scalar reduction operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/detail.hpp"
#include "viennacl/generator/utils.hpp"

#include "viennacl/generator/generate_template_base.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    class scalar_reduction : public template_base{
        typedef template_base base_type;
      public:
        typedef base_type::statements_type statements_type;

        class profile : public template_base::profile{
            friend class scalar_reduction;
            std::size_t lmem_used(std::size_t scalartype_size) const {
              return group_size_*scalartype_size;
            }
          public:
            /** @brief The user constructor */
            profile(unsigned int vectorization, unsigned int group_size, unsigned int num_groups, bool global_decomposition) : template_base::profile(vectorization), group_size_(group_size), num_groups_(num_groups), global_decomposition_(global_decomposition){ }
            void set_local_sizes(std::size_t& s1, std::size_t& s2) const{
              s1 = group_size_;
              s2 = 1;
            }
            void kernel_arguments(std::string & arguments_string) const{
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }
          private:
            unsigned int group_size_;
            unsigned int num_groups_;
            bool global_decomposition_;
        };



      private:

        void core_0(utils::kernel_generation_stream& stream) const {

          std::vector<detail::mapped_scalar_reduction*> exprs;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::iterator iit = it->begin() ; iit != it->end() ; ++iit)
              if(detail::mapped_scalar_reduction * p = dynamic_cast<detail::mapped_scalar_reduction*>(iit->second.get()))
                exprs.push_back(p);
          std::string scalartype = exprs.front()->scalartype();

          for(std::size_t k = 0 ; k < exprs.size() ; ++k)
              stream << scalartype << " sum" << k << " = 0;" << std::endl;



          if(profile_.global_decomposition_){
            stream << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0)){" << std::endl;
          }
          else{
            stream << "unsigned int chunk_size = (N + get_num_groups(0)-1)/get_num_groups(0);" << std::endl;
            stream << "unsigned int chunk_start = get_group_id(0)*chunk_size;" << std::endl;
            stream << "unsigned int chunk_end = min(chunk_start+chunk_size, N);" << std::endl;
            stream << "for(unsigned int i = chunk_start + get_local_id(0) ; i < chunk_end ; i += get_local_size(0)){" << std::endl;
          }
          stream.inc_tab();

          //Fetch vector entry
          std::set<std::string>  fetched;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::reverse_iterator it2 = it->rbegin() ; it2 != it->rend() ; ++it2)
              if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(it2->second.get()))
                p->fetch( "i", fetched, stream);

          //Update sums;
          for(std::size_t k = 0 ; k < exprs.size() ; ++k){
            std::string expr_str;
            detail::traverse(*exprs[k]->lhs().array_, detail::expression_generation_traversal("", expr_str, *exprs[k]->lhs().mapping_), false, exprs[k]->lhs().index_);
            expr_str += "*";
            detail::traverse(*exprs[k]->rhs().array_, detail::expression_generation_traversal("", expr_str, *exprs[k]->rhs().mapping_), false, exprs[k]->rhs().index_);
            stream << scalartype << " sum" << k << " += "  << expr_str << std::endl;
           }

          stream.dec_tab();
          stream << "}" << std::endl;
          //Declare and fill local memory
          for(std::size_t k = 0 ; k < exprs.size() ; ++k){
              stream << "__local " << scalartype << " buf" << k << "[" << profile_.group_size_ << "];" << std::endl;
              stream << "buf" << k << "[lid] = sum" << k << std::endl;
          }
          //Reduce local memory
          for(unsigned int stride = profile_.group_size_/2 ; stride>0 ; stride /=2){
            stream << "barrier(CLK_LOCAL_MEM_FENCE); " << std::endl;
            for(std::size_t k = 0 ; k < exprs.size() ; ++k){
                stream << "if(lid < " << stride << ")" << std::endl;
                stream << "buf" << k << "[lid] += buf" << k << "[lid + " << stride << "];" << std::endl;
            }
          }

          //Fetches to temporary buffer
          stream << "if(lid==0){" << std::endl;
          stream.inc_tab();
          for(std::size_t k = 0 ; k < exprs.size() ; ++k){
              stream << "temp"<< "[get_group_id(0)] = local" << k << "[0];" << std::endl;
          }
          stream.dec_tab();
          stream << "}" << std::endl;
        }


        void core_1(utils::kernel_generation_stream& stream) const {
//          for(std::size_t i = 0 ; i < inner_products_index.size() ; ++i){
//            for(std::size_t j = 0 ; j < inner_products_index[i].size() ; ++j){
//              stream << "__local " << scalartype << " local" << i << j << "[" << profile.group_size_ << "];" << std::endl;
//              stream << "local" << i << j << "[lid] = " << name << "[lid];" << std::endl;
//            }
//          }
          stream << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
//        if(lid==0) " << (*it)->generate(0) << ";" << std::endl;
        }

      public:
        scalar_reduction(template_base::statements_type const & s, profile const & p) : template_base(s, profile_), profile_(p){ }

        void core(std::size_t kernel_id, utils::kernel_generation_stream& stream) const {
          assert(kernel_id<1 && bool("Core not implemented"));
          if(kernel_id==0){
            core_0(stream);
          }
          else{
            core_1(stream);
          }
        }

      private:
        profile profile_;
    };

  }

}

#endif
