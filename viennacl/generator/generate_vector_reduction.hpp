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

          std::size_t lsize1 = profile_.m_;
          std::size_t lsize2 = profile_.k_+1;
          bool is_lhs_transposed = false;

          std::vector<detail::mapped_vector_reduction*> exprs;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::iterator iit = it->begin() ; iit != it->end() ; ++iit)
              if(detail::mapped_vector_reduction * p = dynamic_cast<detail::mapped_vector_reduction*>(iit->second.get()))
                exprs.push_back(p);

          detail::mapped_vector_reduction* first_expr = exprs.front();
          for(std::vector<detail::mapped_vector_reduction*>::iterator it = exprs.begin() ; it != exprs.end() ; ++it){
            stream << "__local " <<  (*it)->scalartype() << " buf" << std::distance(exprs.begin(), it) << '[' << lsize1*lsize2 << ']' << std::endl;
          }

          stream << "unsigned int lid0 = get_local_id(0);" << std::endl;
          stream << "unsigned int lid1 = get_local_id(1);" << std::endl;


          if(is_lhs_transposed)
            stream << "for(unsigned int r = get_global_id(0) ; r < N ; r += get_global_size(0)){" << std::endl;
          else
            stream << "for(unsigned int r = get_global_id(0) ; r < M ; r += get_global_size(0)){" << std::endl;
          stream.inc_tab();

          stream << first_expr->scalartype() << " sum = 0;" << std::endl;
          if(is_lhs_transposed)
            stream << "for(unsigned int c = get_local_id(1) ; c < M ; c += get_local_size(1)){" << std::endl;
          else
            stream << "for(unsigned int c = get_local_id(1) ; c < N ; c += get_local_size(1)){" << std::endl;
          stream.inc_tab();

          std::string expr_str;
          detail::traverse(*first_expr->lhs().array_, detail::expression_generation_traversal("r*N+c", expr_str, *first_expr->lhs().mapping_), false, first_expr->lhs().index_);
          expr_str += "*";
          detail::traverse(*first_expr->rhs().array_, detail::expression_generation_traversal("c", expr_str, *first_expr->rhs().mapping_), false, first_expr->rhs().index_);

          stream << "sum = sum" << "+" << '(' << expr_str << ");" << std::endl;
          stream.dec_tab();
          stream << "}" << std::endl;

          stream << "buf[lid0*" << lsize2 << "+ lid1] = sum" << std::endl;

          for(unsigned int stride = profile_.k_/2 ; stride>1 ; stride /=2){
            stream << "barrier(CLK_LOCAL_MEM_FENCE); " << std::endl;
            stream <<  "if(lid1 < " << stride << ")" ;
            stream << "{" << std::endl;
            stream.inc_tab();

            //Compute final expression in final reduction
            for(std::size_t i = 0 ; i < exprs.size() ; ++i)
              stream << "buf" << i << "[lid0*" << lsize2 << "+ lid1] += buf" << i << "[lid0*" << lsize2 << "+ " << stride << "];" << std::endl;

            stream.dec_tab();
            stream << "}" << std::endl;
          }


          stream << "barrier(CLK_LOCAL_MEM_FENCE); " << std::endl;
          stream <<  "if(lid1 == 0)" ;
          stream << "{" << std::endl;
          stream.inc_tab();
          for(std::size_t i = 0 ; i < exprs.size() ; ++i){
            stream << "buf" << i << "[lid0*" << lsize2 << "] += buf" << i << "[lid0*" << lsize2 << "+ 1];" << std::endl;
            exprs[i]->access_name("buf"+utils::to_string(i)+"[lid0*"+utils::to_string(lsize2)+"]");
          }
          std::size_t i = 0;
          for(statements_type::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it){
            std::string str;
            detail::traverse(it->array(), detail::expression_generation_traversal("r", str, mapping_[i++]), false);
            stream << str << std::endl;
          }
          stream.dec_tab();
          stream << "}" << std::endl;


        }

      private:
        profile profile_;
    };
  }
}

#endif
