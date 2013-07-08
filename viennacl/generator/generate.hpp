#ifndef VIENNACL_GENERATOR_GENERATE_HPP
#define VIENNACL_GENERATOR_GENERATE_HPP

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


/** @file viennacl/generator/generate.hpp
    @brief User interface
*/

#include <vector>
#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/detail.hpp"
#include "viennacl/generator/generate_saxpy_vector.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{

      class name_generation_traversal : public traversal_functor{
          std::string & str_;
        public:
          name_generation_traversal(std::string & str) : str_(str) { }
          void call_on_leaf(std::size_t, statement_node_type_family, statement_node_type type, lhs_rhs_element) const { str_ += detail::generate(type); }
          void call_on_op(operation_node_type_family, operation_node_type type) const { str_ += detail::generate(type); }
          void call_before_expansion() const { str_ += '('; }
          void call_after_expansion() const { str_ += ')'; }
      };

      class prototype_generation_traversal : public traversal_functor{

          mutable unsigned int current_arg_;
          std::map<cl_mem, std::size_t> & memory_;
          mapping_type & mapping_;
          std::string & str_;


          void prototype_value_generation(statement_node_type_family, statement_node_type type, lhs_rhs_element, symbolic_container & sym) const{
            sym.scalartype_ = detail::generate_scalartype(type);
            sym.name_ = "arg" + utils::to_string(current_arg_++);
            str_ += sym.scalartype_ + ' '  + sym.name_ + ",";
          }

          void prototype_pointer_generation(statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element, symbolic_container & sym) const {
            sym.scalartype_ = detail::generate_scalartype(type);
            if(memory_.insert(std::make_pair(detail::get_handle(type, element), current_arg_)).second){
              sym.name_ =  "arg" + utils::to_string(current_arg_++);
              str_ += "__global " +  sym.scalartype_ + "* " + sym.name_ + ",";
            }
            sym.name_ = "arg" + utils::to_string(memory_.at(detail::get_handle(type, element)));
          }

        public:
          prototype_generation_traversal(std::map<cl_mem, std::size_t> & memory, mapping_type & mapping, std::string & str) : current_arg_(0), memory_(memory), mapping_(mapping), str_(str) { }

          void call_on_leaf(std::size_t index, statement_node_type_family type_family, statement_node_type type, lhs_rhs_element element) const {
            if(type_family==HOST_SCALAR_TYPE_FAMILY)
              prototype_value_generation(type_family,type,element, mapping_[index]);
            else
              prototype_pointer_generation(type_family,type,element, mapping_[index]);
          }
      };
    }

    using namespace viennacl::scheduler;

    std::string make_program_name(std::vector<scheduler::statement> const & statements){
      std::string res;
      for(std::vector<scheduler::statement>::const_iterator it = statements.begin() ; it != statements.end() ; ++it)
        detail::traverse(it->array(), 0, detail::name_generation_traversal(res));
      return res;
    }

    std::string make_program_string(std::vector<scheduler::statement> const & statements){
      std::map<cl_mem, std::size_t> memory;
      std::size_t size = statements.size();
      std::vector<detail::mapping_type> mapping(size);

      utils::kernel_generation_stream kss;

      //Headers generation
      kss << "#if defined(cl_khr_fp64)\n";
      kss <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
      kss <<  "#elif defined(cl_amd_fp64)\n";
      kss <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
      kss <<  "#endif\n";

      kss << std::endl;

      //Prototype generation
      kss << "__kernel void kernel_0(" << std::endl;
      std::string prototype;
      for(std::size_t i = 0 ; i < size ; ++i)
        detail::traverse(statements[i].array(), 0, detail::prototype_generation_traversal(memory, mapping[i], prototype));
      prototype.erase(prototype.size()-1); //Last comma pruned
      kss << prototype << std::endl;
      kss << ")" << std::endl;
      kss << "{" << std::endl;
      kss.inc_tab();

      //body generation
      {
        generate_saxpy_vector(saxpy_vector_profile(1,4,128), kss, statements, mapping);
      }

      kss.dec_tab();
      kss << "}" << std::endl;
      return kss.str();
    }

  }

}
#endif
