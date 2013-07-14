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

    class saxpy{
      public:
        typedef std::vector<scheduler::statement> statements_type;

        class profile : public profile_base{
          public:
            profile(unsigned int vectorization, std::size_t group_size, std::size_t num_groups, bool global_decomposition) : profile_base(vectorization)
                                                                                                    , group_size_(group_size)
                                                                                                    , num_groups_(num_groups)
                                                                                                    , global_decomposition_(global_decomposition){ }

            void set_local_sizes(std::size_t & x, std::size_t & y) const{
              x = group_size_;
              y = 1;
            }
            static void kernel_arguments(std::string & arguments_string){
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }
          private:
            std::size_t group_size_;
            std::size_t num_groups_;
            bool global_decomposition_;
        };

      private:
        profile profile_;
        statements_type const & statements_;
        mutable std::vector<detail::mapping_type> mapping_;

        void fetch(detail::mapped_handle * c, std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) const {
          std::string new_access_name = c->name() + "_private";
          if(fetched.find(c->name())==fetched.end()){
            stream << c->scalartype() << " " << new_access_name << " = ";
            c->generate(index, stream);
            stream << ';' << std::endl;
            fetched.insert(c->name());
          }
          c->access_name(new_access_name);
        }

        void write_back(detail::mapped_handle * c, std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) const {
          std::string old_access_name = c->access_name();
          c->access_name("");
          if(fetched.find(c->name())!=fetched.end()){
            c->generate(index, stream);
            stream << " = " << old_access_name << ';' << std::endl;
            fetched.erase(c->name());
          }
        }
      public:
        saxpy(profile const & p, statements_type const & s) : profile_(p), statements_(s), mapping_(s.size()){ }

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

        void core(utils::kernel_generation_stream& stream) const {
          stream << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0))" << std::endl;
          stream << "{" << std::endl;
          stream.inc_tab();

    //      if(first_matrix->is_rowmajor()){
    //        kss << "unsigned int r = get_global_id(0)/" << first_matrix->internal_size2() << ";" << std::endl;
    //        kss << "unsigned int c = get_global_id(0)%" << first_matrix->internal_size2() << ";" << std::endl;
    //      }
    //      else{
    //        kss << "unsigned int r = get_global_id(0)%" << first_matrix->internal_size1() << ";" << std::endl;
    //        kss << "unsigned int c = get_global_id(0)/" << first_matrix->internal_size1() << ";" << std::endl;
    //      }

          //Fetches entries to registers
          std::set<std::string>  fetched;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::reverse_iterator it2 = it->rbegin() ; it2 != it->rend() ; ++it2)
              if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(it2->second.get()))
                fetch(p, "i", fetched, stream);

          std::size_t i = 0;
          for(statements_type::const_iterator it = statements_.begin() ; it != statements_.end() ; ++it){
              detail::traverse(it->array(), detail::expression_generation_traversal("i", stream,mapping_[i++]), true);
              stream << ";" << std::endl;
          }

          //Writes back
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(it->at(std::make_pair(0,detail::LHS_LEAF_TYPE)).get()))
              write_back(p, "i", fetched, stream);

          stream.dec_tab();
          stream << "}" << std::endl;
        }
    };

  }

}

#endif