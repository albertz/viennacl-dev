#ifndef VIENNACL_GENERATOR_GENERATE_SAXPY_HPP
#define VIENNACL_GENERATOR_GENERATE_SAXPY_HPP

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

#include "viennacl/generator/generate_utils.hpp"
#include "viennacl/generator/utils.hpp"

#include "viennacl/generator/generate_template_base.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    class vector_saxpy : public template_base{
      public:
        class profile : public template_base::profile{
          public:
            profile(unsigned int v, std::size_t gs, std::size_t ng, bool d) : template_base::profile(v), group_size_(gs), num_groups_(ng), global_decomposition_(d){ }
            void set_local_sizes(std::size_t & size1, std::size_t & size2) const{
              size1 = group_size_;
              size2 = 1;
            }
            virtual void enqueue_kernel_arguments(statements_type  const & statements, viennacl::ocl::kernel & k, unsigned int & n_arg, unsigned int kernel_id) const{
              scheduler::statement_node first_node = statements.front().array()[0];
              k.arg(n_arg++, cl_uint(utils::call_on_vector(first_node.lhs_type_, first_node.lhs_, utils::size_fun())));
            }
            void kernel_arguments(statements_type  const & statements, std::string & arguments_string) const{
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }
          private:
            std::size_t group_size_;
            std::size_t num_groups_;
            bool global_decomposition_;
        };

      public:
        vector_saxpy(template_base::statements_type const & s, profile const & p) : template_base(s, 1, profile_), profile_(p){ }

        void core(std::size_t kernel_id, utils::kernel_generation_stream& stream) const {
          stream << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0))" << std::endl;
          stream << "{" << std::endl;
          stream.inc_tab();

          //Fetches entries to registers
          std::set<std::string>  fetched;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::reverse_iterator it2 = it->rbegin() ; it2 != it->rend() ; ++it2)
              if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(it2->second.get()))
                p->fetch( "i", fetched, stream);

          std::size_t i = 0;
          for(statements_type::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it){
              std::string str;
              detail::traverse(it->array(), detail::expression_generation_traversal("", str, mapping_[i++]), false);
              stream << str << std::endl;
          }

          //Writes back
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(it->at(std::make_pair(0,detail::LHS_NODE_TYPE)).get()))
              p->write_back( "i", fetched, stream);

          stream.dec_tab();
          stream << "}" << std::endl;
        }

      private:
        profile profile_;
    };



    class matrix_saxpy : public template_base{
      public:
        class profile : public template_base::profile{
          public:
            profile(unsigned int v, std::size_t gs1, std::size_t gs2, std::size_t ng1, std::size_t ng2, bool d) : template_base::profile(v), group_size1_(gs1), group_size2_(gs2), num_groups1_(ng1), num_groups2_(ng2), global_decomposition_(d){ }
            void set_local_sizes(std::size_t & x, std::size_t & y) const{
              x = group_size1_;
              y = group_size2_;
            }
            void enqueue_kernel_arguments(statements_type  const & statements, viennacl::ocl::kernel & k, unsigned int & n_arg, unsigned int kernel_id) const{
              scheduler::statement_node first_node = statements.front().array()[0];
              k.arg(n_arg++, cl_uint(utils::call_on_matrix(first_node.lhs_type_, first_node.lhs_, utils::size1_fun())));
              k.arg(n_arg++, cl_uint(utils::call_on_matrix(first_node.lhs_type_, first_node.lhs_, utils::size2_fun())));
            }
            void kernel_arguments(statements_type  const & statements, std::string & arguments_string) const{
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "M");
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }
          private:
            std::size_t group_size1_;
            std::size_t group_size2_;

            std::size_t num_groups1_;
            std::size_t num_groups2_;

            bool global_decomposition_;
        };

      private:
        std::string get_offset(detail::mapped_matrix * p) const {
          if(p->is_row_major())
            return "i*N + j";
          else
            return "i + j*M";
        }

      public:
        matrix_saxpy(template_base::statements_type const & s, profile const & p) : template_base(s, 1, profile_), profile_(p){ }

        void core(std::size_t kernel_id, utils::kernel_generation_stream& stream) const {
          stream << "for(unsigned int i = get_global_id(0) ; i < M ; i += get_global_size(0))" << std::endl;
          stream << "{" << std::endl;
          stream.inc_tab();
          stream << "for(unsigned int j = get_global_id(1) ; j < N ; j += get_global_size(1))" << std::endl;
          stream << "{" << std::endl;
          stream.inc_tab();


          //Fetches entries to registers
          std::set<std::string>  fetched;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::reverse_iterator it2 = it->rbegin() ; it2 != it->rend() ; ++it2)
              if(detail::mapped_matrix * p = dynamic_cast<detail::mapped_matrix *>(it2->second.get()))
                p->fetch(get_offset(p),fetched, stream);

          std::size_t i = 0;
          for(statements_type::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it){
            std::string str;
            detail::traverse(it->array(), detail::expression_generation_traversal("", str, mapping_[i++]), false);
            stream << str << std::endl;
          }

          //Writes back
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            if(detail::mapped_matrix * p = dynamic_cast<detail::mapped_matrix *>(it->at(std::make_pair(0,detail::LHS_NODE_TYPE)).get()))
              p->write_back(get_offset(p), fetched, stream);

          stream.dec_tab();
          stream << "}" << std::endl;
          stream.dec_tab();
          stream << "}" << std::endl;
        }

      private:
        profile profile_;
    };
  }

}

#endif
