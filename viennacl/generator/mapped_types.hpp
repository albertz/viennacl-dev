#ifndef VIENNACL_GENERATOR_MAPPED_TYPE_HPP
#define VIENNACL_GENERATOR_MAPPED_TYPE_HPP

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


/** @file viennacl/generator/mapped_types.hpp
    @brief Internal upper layer to the scheduler
*/

#include <string>

#include "viennacl/scheduler/forwards.h"
#include "viennacl/generator/forwards.h"
#include "viennacl/generator/utils.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{


      /** @brief Base class for mapping viennacl datastructure to generator-friendly structures
       */
      class mapped_container{
        protected:
          struct node_info{
              node_info() : mapping_(NULL), array_(NULL) { }
              mapping_type const * mapping_;
              statement::container_type const * array_;
              index_info index_;
          };
          virtual std::string generate_default(std::string const & index) const = 0;
        public:
          mapped_container(std::string const & scalartype) : scalartype_(scalartype){          }
          std::string const & scalartype() const { return scalartype_; }
          void access_name(std::string const & str) { access_name_ = str; }
          std::string const & access_name() const { return access_name_; }
          virtual std::string generate(std::string const & index) const{
            if(!access_name_.empty())
              return access_name_;
            else
              return generate_default(index);
          }
          virtual std::string & append_kernel_arguments(std::string & str){ }
          virtual ~mapped_container(){ }
        protected:
          std::string access_name_;
          std::string scalartype_;
      };

      class mapped_binary_leaf : public mapped_container{
        public:
          mapped_binary_leaf(std::string const & scalartype) : mapped_container(scalartype){ }
          node_info lhs() const { return lhs_; }
          node_info rhs() const { return rhs_; }
          std::string generate_default(std::string const & index) const { return "";}
        protected:
          node_info lhs_;
          node_info rhs_;
      };

      class mapped_matrix_product : public mapped_binary_leaf{
          friend class map_generate_prototype;
        public:
          mapped_matrix_product(std::string const & scalartype) : mapped_binary_leaf(scalartype){ }
      };

      class mapped_reduction : public mapped_binary_leaf{
        public:
          mapped_reduction(std::string const & scalartype) : mapped_binary_leaf(scalartype){ }
          operation_node_type reduction_type() const { return reduction_type_; }
        private:
          operation_node_type reduction_type_;
      };

      class mapped_scalar_reduction : public mapped_reduction{
          friend class map_generate_prototype;
        public:
          mapped_scalar_reduction(std::string const & scalartype) : mapped_reduction(scalartype){ }
      };

      class mapped_vector_reduction : public mapped_reduction{
          friend class map_generate_prototype;
        public:
          mapped_vector_reduction(std::string const & scalartype) : mapped_reduction(scalartype){ }
      };

      /** @brief Mapping of a host scalar to a generator class */
      class mapped_host_scalar : public mapped_container{
          friend class map_generate_prototype;
          std::string generate_default(std::string const & index) const{ return name_;  }
        public:
          mapped_host_scalar(std::string const & scalartype) : mapped_container(scalartype){ }
          std::string const & name() { return name_; }
        private:
          std::string name_;
      };

      /** @brief Base class for datastructures passed by pointer */
      class mapped_handle : public mapped_container{
          virtual std::string offset(std::string const & index) const = 0;
          std::string generate_default(std::string const & index) const{ return name_  + '[' + offset(index) + ']'; }
        public:
          mapped_handle(std::string const & scalartype) : mapped_container(scalartype){ }

          std::string const & name() const { return name_; }

          void fetch(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) {
            std::string new_access_name = name_ + "_private";
            if(fetched.find(name_)==fetched.end()){
              stream << scalartype_ << " " << new_access_name << " = " << generate_default(index) << ';' << std::endl;
              fetched.insert(name_);
            }
            access_name_ = new_access_name;
          }

          void write_back(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) {
            std::string old_access_name = access_name_ ;
            access_name_ = "";
            if(fetched.find(name_)!=fetched.end()){
              stream << generate_default(index) << " = " << old_access_name << ';' << std::endl;
              fetched.erase(name_);
            }
          }
        protected:
          std::string name_;
      };

      /** @brief Base class for scalar */
      class mapped_scalar : public mapped_handle{
          friend class map_generate_prototype;
        private:
          std::string offset(std::string const & index)  const { return "0"; }
        public:
          mapped_scalar(std::string const & scalartype) : mapped_handle(scalartype){ }
      };


      /** @brief Mapping of a vector to a generator class */
      class mapped_vector : public mapped_handle{
          friend class map_generate_prototype;
          std::string offset(std::string const & index) const {
            if(access_node_.array_){
              std::string str;
              detail::traverse(*access_node_.array_, detail::expression_generation_traversal(index, str, *access_node_.mapping_), true, access_node_.index_);
              return str;
            }
            else
              return index;
          }
        public:
          mapped_vector(std::string const & scalartype) : mapped_handle(scalartype){ }
        private:
          node_info access_node_;
          std::string start_name_;
          std::string stride_name_;
          std::string shift_name_;
      };

      /** @brief Mapping of a matrix to a generator class */
      class mapped_matrix : public mapped_handle{
          friend class map_generate_prototype;
          std::string offset(std::string const & index) const {
            return index;
          }
        public:
          bool is_row_major() const { return is_row_major_; }
          bool is_transposed() const { return is_transposed_; }
          std::string const & size1() const { return size1_; }
          std::string const & size2() const { return size2_; }
          std::string offset(std::string const & i, std::string const & j) const{
            if(is_row_major_)
              if(j=="0")
                return '(' + i + ')' + '*' + size2_;
              else
                return '(' + i + ')' + '*' + size2_ + "+ (" + j + ')';
            else
              if(i=="0")
                return  "(" + j + ')' + '*' + size1_;
              else
                return  '(' + i + ')' + "+ (" + j + ')' + '*' + size1_;
          }
          mapped_matrix(std::string const & scalartype) : mapped_handle(scalartype){ }
        private:
          std::string size1_;
          std::string size2_;

          std::string start1_name_;
          std::string stride1_name_;
          std::string shift1_name_;
          std::string start2_name_;
          std::string stride2_name_;
          std::string shift2_name_;
          bool is_row_major_;
          bool is_transposed_;
      };

      /** @brief Mapping of a symbolic vector to a generator class */
      class mapped_symbolic_vector : public mapped_container{
          friend class map_generate_prototype;
          std::string value_name_;
          std::string index_name_;
          bool is_value_static_;
        public:
          mapped_symbolic_vector(std::string const & scalartype) : mapped_container(scalartype){ }
          std::string generate_default(std::string const & index) const{
            return value_name_;
          }
      };

      /** @brief Mapping of a symbolic matrix to a generator class */
      class mapped_symbolic_matrix : public mapped_container{
          friend class map_generate_prototype;
          std::string value_name_;
          bool is_diag_;
        public:
          mapped_symbolic_matrix(std::string const & scalartype) : mapped_container(scalartype){ }
          std::string generate_default(std::string const & index) const{
            return value_name_;
          }
      };

      std::string generate(std::string const & index, mapped_container const & s){
        return s.generate(index);
      }

      void fetch(std::string const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream, mapped_container & s){
        if(mapped_handle * p = dynamic_cast<mapped_handle  *>(&s))
          p->fetch(index, fetched, stream);
      }

    }

  }

}
#endif
