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

#include "viennacl/generator/template_base.hpp"

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

            void kernel_arguments(std::string & arguments_string) const{
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }

          private:
            unsigned int group_size_;
            unsigned int num_groups_;
            bool global_decomposition_;
        };

        void scalar_reduction_1(utils::kernel_generation_stream& stream) const {
          std::vector< std::vector<std::size_t> > inner_products_index;
          std::string scalartype = mapping_[0][std::make_pair(0,detail::LHS_LEAF_TYPE)];

          for(statements_type::iterator it = statements_.begin() ; it != statements_.end() ; ++it){
            std::vector<std::size_t> indexes;
            for(scheduler::statement::container_type::iterator iit = it->array().begin() ; iit != it->array.end() ; ++iit)
              if(iit->op_family_ = OPERATION_REDUCE_TYPE_FAMILY)
                indexes.push_back(std::distance(it->array().begin(), iit));
            inner_products_index.push_back(indexes);
          }


          for(std::size_t i = 0 ; i != inner_products_index.size() ; ++i)
            for(std::size_t j = 0 ; j != inner_products_index[i].size() ; ++j)
              kss << scalartype << " sum" << i << j << " = 0;" ;



          if(profile_.global_decomposition_){
            kss << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0)){" << std::endl;
          }
          else{
            kss << "unsigned int chunk_size = (N + get_num_groups(0)-1)/get_num_groups(0);" << std::endl;
            kss << "unsigned int chunk_start = get_group_id(0)*chunk_size;" << std::endl;
            kss << "unsigned int chunk_end = min(chunk_start+chunk_size, " << size << ");" << std::endl;
            kss << "for(unsigned int i = chunk_start + get_local_id(0) ; i < chunk_end ; i += get_local_size(0)){" << std::endl;
          }
          kss.inc_tab();

          //Fetch
          std::set<std::string>  fetched;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::reverse_iterator it2 = it->rbegin() ; it2 != it->rend() ; ++it2)
              if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(it2->second.get()))
                p->fetch( "i", fetched, stream);


          //Update sums;
          for(std::size_t i = 0 ; i < inner_products_index.size() ; ++i)
            for(std::size_t j = 0 ; j < inner_products_index[i].size() ; ++j)
              for(unsigned int v = 0 ; v < profile_.vectorization_ ; ++v)
                kss << scalartype << " sum" << i << j << " += "  << detail::traverse(statements_[i].array(), detail::expression_generation_traversal("", stream, mapping_[j], true));

          kss.dec_tab();
          kss << "}" << std::endl;

          //Declare and fill local memory
          for(std::size_t i = 0 ; i < inner_products_index.size() ; ++i){
            for(std::size_t j = 0 ; j < inner_products_index[i].size() ; ++j){
              kss << "__local " << scalartype << " local" << i << j << "[" << profile.group_size_ << "];" << std::endl;
              kss << "local" << i << j << "[lid] = sum" << i << j << std::endl;
            }
          }

          //Reduce local memory
          for(unsigned int stride = profile.group_size_/2 ; stride>0 ; stride /=2){
            kss << "barrier(CLK_LOCAL_MEM_FENCE); " << std::endl;
            for(std::size_t i = 0 ; i < inner_products_index.size() ; ++i){
              for(std::size_t j = 0 ; j < inner_products_index[i].size() ; ++j){
                kss << "local" << i << j << "[lid] = (lid < " << stride << ")? local" << i << j << "[lid + " << stride << "]:0;" << std::endl;
              }
            }
          }

          for(std::size_t i = 0 ; i < inner_products_index.size() ; ++i){
            for(std::size_t j = 0 ; j < inner_products_index[i].size() ; ++j){
              kss << "if(lid==0) local" << i << j << "[get_group_id(0)] = local" << i << j << "[0];" << std::endl;
            }
          }
        }

        void scalar_reduction_2(utils::kernel_generation_stream& stream) const {
          for(std::size_t i = 0 ; i < inner_products_index.size() ; ++i){
            for(std::size_t j = 0 ; j < inner_products_index[i].size() ; ++j){
              kss << "__local " << scalartype << " local" << i << j << "[" << profile.group_size_ << "];" << std::endl;
              kss << "local" << i << j << "[lid] = " << name << "[lid];" << std::endl;
            }
          }
          kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
//        if(lid==0) " << (*it)->generate(0) << ";" << std::endl;
        }

      public:
        scalar_reduction( const & s, profile const & p) : template_base(s, profile_), profile_(p){ }

        void core(utils::kernel_generation_stream& stream) const{

        }
      private:
        profile profile_;
    };

}

#endif
