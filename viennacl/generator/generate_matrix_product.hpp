#ifndef VIENNACL_GENERATOR_GENERATE_MATRIX_PRODUCT_HPP
#define VIENNACL_GENERATOR_GENERATE_MATRIX_PRODUCT_HPP

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


/** @file viennacl/generator/templates/matrix_product.hpp
 *
 * Kernel template for the vector reduction operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/generate_template_base.hpp"
#include "viennacl/generator/mapped_types.hpp"
#include "viennacl/generator/utils.hpp"


#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    class matrix_product : public template_base{
        typedef template_base base_type;
      public:
        typedef base_type::statements_type statements_type;

        class profile : public template_base::profile{
            friend class matrix_product;
            std::size_t lmem_used(std::size_t scalartype_size) const {
              std::size_t lmem_used = 0;
              if(use_LHS_shared_) lmem_used += (ml_ + 1) * (kl_ + 1) * scalartype_size;
              if(use_RHS_shared_) lmem_used += (kl_ + 1) * (nl_ + 1) * scalartype_size;
              return lmem_used;
            }

          public:
            /** @brief The user constructor */
            profile(unsigned int vectorization, unsigned int ml, unsigned int kl, unsigned int nl
                    , unsigned int ms, unsigned int ks, unsigned int ns
                    , bool use_LHS_shared, bool use_RHS_shared
                    , unsigned int unroll) : template_base::profile(vectorization){
              ml_= ml; kl_=kl ; nl_=nl;
              ms_ = ms; ks_=ks; ns_=ns;
              use_LHS_shared_ = use_LHS_shared ; use_RHS_shared_ = use_RHS_shared;
              vectorization_ = vectorization;
              unroll_ = unroll;
            }

            void set_local_sizes(std::size_t& s1, std::size_t& s2) const{
              s1 = ml_/ms_;
              s2 = nl_/ns_;
            }

            static std::string size1() { return "M";  }
            static std::string size2() { return "K"; }
            static std::string size3() { return "N"; }

            void kernel_arguments(statements_type  const & statements, std::string & arguments_string) const{
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "M");
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "K");
              arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
            }

          private:
            unsigned int ml_;
            unsigned int kl_;
            unsigned int nl_;

            unsigned int ms_;
            unsigned int ks_;
            unsigned int ns_;

            bool use_LHS_shared_;
            bool use_RHS_shared_;

            unsigned int unroll_;
        };

        static void transform_block(bool is_rowmajor, std::size_t vectorization
                                    , bool is_transposed, bool store_shared
                                    , unsigned int & large_block_1, unsigned int & large_block_2
                                    , unsigned int & small_block_1, unsigned int & small_block_2){
          if(is_rowmajor){
            if(is_transposed) large_block_1 /= vectorization;
            else large_block_2/=vectorization;
            if(!store_shared){
              if(is_transposed) small_block_1/=vectorization;
              else small_block_2/=vectorization;
            }
          }
          else{
            if(is_transposed) large_block_2 /= vectorization;
            else large_block_1/=vectorization;
            if(!store_shared){
              if(is_transposed)  small_block_2/=vectorization;
              else    small_block_1/=vectorization;
            }
          }
        }

        std::string helper_variable(utils::kernel_generation_stream & stream
                                    , bool store_in_register
                                    , std::string const & type
                                    , std::string const & name
                                    , std::string const & expr){
          if(!store_in_register)
            return expr;
          stream << type << " " << name << " = " << expr << ";" << std::endl;
          return name;
        }

        void init_rhs_global_ptr(utils::kernel_generation_stream &  stream, detail::mapped_matrix const & mat, std::string const & offset, unsigned int ks_rhs,unsigned int ns_rhs, unsigned int nl_rhs) {
          if(mat.is_row_major())
            for(unsigned int k = 0 ; k < ks_rhs ; ++k){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(k);
              stream << "__global " << mat.scalartype() << profile_.vectorization_ << " * " << ptr_name << " = " << mat.name() << " + " ;
              if(mat.is_transposed()){
                std::string i = utils::to_string(k) + " + " + offset + " +  get_group_id(1)*" + utils::to_string(nl_rhs);
                stream << mat.offset(i, "0");
              }
              else{
                std::string i = utils::to_string(k);
                std::string j = offset + " +  get_group_id(1)*" + utils::to_string(nl_rhs);
                stream << mat.offset(i, j);
              }
              stream << ";" << std::endl;
              //            mat->private_value(k,"*" + ptr_name);
            }
          else
            for(unsigned int n = 0 ; n < ns_rhs ; ++n){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(n);
              stream << "__global " << mat.scalartype() << profile_.vectorization_ << " * " << ptr_name << " = " << mat.name() << " +  " ;
              if(mat.is_transposed()){
                std::string i = offset + " +  get_group_id(1)*" + utils::to_string(nl_rhs);
                std::string j = utils::to_string(n);
                stream << mat.offset(i, j);
              }
              else{
                std::string j = offset + " +  get_group_id(1)*" + utils::to_string(nl_rhs) + " + " + utils::to_string(n);
                stream << mat.offset("0",j);
              }
              stream << ";" << std::endl;
              //            mat->private_value(n,"*" + ptr_name);
            }
        }

        void update_rhs_global_ptr(utils::kernel_generation_stream & stream, detail::mapped_matrix const & mat, std::string const & offset, unsigned int ks, unsigned int ns_rhs, unsigned int ks_rhs){
          if(mat.is_row_major() && !mat.is_transposed())
            for(unsigned int k=0 ; k<ks ; ++k)
              stream << mat.name() << "_ptr_" << k << " += " << ks_rhs << "*" << profile::size3() << " - " << ns_rhs << ";" << std::endl;
          else if(mat.is_transposed() && !mat.is_row_major())
            for(unsigned int n=0 ; n<ns_rhs ; ++n)
              stream << mat.name() << "_ptr_" << n << " += " << ns_rhs << "*" << profile::size2() << " - " << ks_rhs << ";" << std::endl;
        }


        void init_lhs_global_ptr(utils::kernel_generation_stream & stream, detail::mapped_matrix const & mat, std::string const & offset, unsigned int ms_lhs, unsigned int ks_lhs, unsigned int ml_lhs) {
          if(mat.is_row_major()){
            for(unsigned int m=0; m<ms_lhs; ++m){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(m);
              stream << "__global " << mat.scalartype() << profile_.vectorization_ << " * " << ptr_name << " = " << mat.name() << " + ";
              if(mat.is_transposed()){
                std::string i = utils::to_string(m);
                std::string j = "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset;
                stream << mat.offset(i,j);
              }
              else{
                std::string i = "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset + "+" + utils::to_string(m);
                stream << mat.offset(i,"0");
              }
              stream << ";" << std::endl;
              //            mat.private_value(m,"*" + ptr_name);
            }
          }
          else{
            for(unsigned int k=0; k<ks_lhs; ++k){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(k);
              stream << "__global " << mat.scalartype() << profile_.vectorization_ << " * " << ptr_name << " = " << mat.name() << " + " ;
              if(mat.is_transposed()){
                std::string j = utils::to_string(k) + "+" + "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset ;
                stream << mat.offset("0",j);
              }
              else{
                std::string i = "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset;
                std::string j = utils::to_string(k);
                stream << mat.offset(i,j);
              }
              stream << ";" << std::endl;
              //            mat.private_value(k,"*" + ptr_name);
            }
          }
        }

        void update_lhs_global_ptr(utils::kernel_generation_stream & stream, detail::mapped_matrix const & mat, unsigned int ks, unsigned int ms_lhs, unsigned int ks_lhs){
          if(mat.is_transposed() && mat.is_row_major())
            for(unsigned int m=0 ; m<ms_lhs ; ++m)
              stream << mat.name() << "_ptr_" << m << " += " << ks << "*" << profile::size2() << " - " <<  ks_lhs << ";" << std::endl;
          else if(!mat.is_transposed() && !mat.is_row_major())
            for(unsigned int k=0 ; k<ks_lhs ; ++k)
              stream << mat.name() << "_ptr_" << k << " += " << ks_lhs << "*" << profile::size1() << " - " << ms_lhs << ";" << std::endl;
        }



        void fetch_to_local_mem(utils::kernel_generation_stream & kss,
                                std::string const & lmem_name,
                                std::size_t lmem_size2,
                                std::string const & offset,
                                unsigned int bound1,
                                unsigned int bound2,
                                detail::mapped_matrix const & mat){
          kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          kss << "for(unsigned int i = get_local_id(0)" << " ; i < " << bound1 << "; i+= get_local_size(0)){" << std::endl;
          kss.inc_tab();
          kss << "for(unsigned int j = get_local_id(1)" << " ; j < " << bound2 << "; j+= get_local_size(1)){" << std::endl;
          kss.inc_tab();
          if(mat.is_row_major()){
            std::string i = "[" + offset + " + j  + " + mat.size2() + "*i]";
            kss << mat.scalartype() << profile_.vectorization_ << " val = " << mat.generate(i) << ";" << std::endl;
            kss << "__local " << mat.scalartype() << "* ptr = " << lmem_name << " + i*" << lmem_size2 << "+j*" << profile_.vectorization_<<";" <<std::endl;
            for(unsigned int a = 0 ; a < profile_.vectorization_ ; ++a){
              if(profile_.vectorization_>1)
                kss << "*ptr++ =  val.s" << a << ";" << std::endl;
              else
                kss << "*ptr++ =  val;" << std::endl;
            }
          }
          else{
            std::string i = mat.name() + "[" + offset + "+ j*" + mat.size1() + " + i]";
            kss << mat.scalartype() << profile_.vectorization_ << " val = " << mat.generate(i) << ";" << std::endl;
            kss << "__local " << mat.scalartype() << "* ptr = " << lmem_name << " + i*" << profile_.vectorization_ * lmem_size2 << "+ j;" <<std::endl;
            for(unsigned int a = 0 ; a < profile_.vectorization_ ; ++a){
              if(profile_.vectorization_>1)
                kss << "*ptr =  val.s" << a << ";" << std::endl;
              else
                kss << "*ptr =  val;" << std::endl;
              kss << "ptr += " << lmem_size2 << ";" << std::endl;
            }
          }

          kss.dec_tab();
          kss << "}" << std::endl;
          kss.dec_tab();
          kss << "}" << std::endl;
          kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        }



      public:
        matrix_product(template_base::statements_type const & s, profile const & p) : template_base(s, profile_), profile_(p){ }

        void core(utils::kernel_generation_stream& stream) const{

          bool use_LHS_shared = profile_.use_LHS_shared_;
          bool use_RHS_shared = profile_.use_RHS_shared_;
          unsigned int vectorization = profile_.vectorization_;
          unsigned int kl = profile_.kl_;
          unsigned int ks = profile_.ks_;
          unsigned int ml = profile_.ml_;
          unsigned int ms = profile_.ms_;
          unsigned int nl = profile_.nl_;
          unsigned int ns = profile_.ns_;
          unsigned int unroll = profile_.unroll_;

          detail::mapped_matrix const * assigned = static_cast<detail::mapped_matrix const *>(mapping_[0].at(std::make_pair(0,detail::LHS_NODE_TYPE)).get());
          std::vector<detail::mapped_matrix_product*> exprs;
          for(std::vector<detail::mapping_type>::iterator it = mapping_.begin() ; it != mapping_.end() ; ++it)
            for(detail::mapping_type::iterator iit = it->begin() ; iit != it->end() ; ++iit)
              if(detail::mapped_matrix_product * p = dynamic_cast<detail::mapped_matrix_product*>(iit->second.get()))
                exprs.push_back(p);
          detail::mapped_matrix const * lhs = static_cast<detail::mapped_matrix const *>(exprs.front()->lhs().mapping_->at(exprs.front()->lhs().index_).get());
          detail::mapped_matrix const * rhs = static_cast<detail::mapped_matrix const *>(exprs.front()->rhs().mapping_->at(exprs.front()->rhs().index_).get());
        }

      private:
        profile profile_;
    };

  }

}

#endif
