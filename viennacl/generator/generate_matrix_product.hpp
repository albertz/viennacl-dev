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

            static std::string size1() { return "M";  }

            static std::string size2() { return "K"; }

            static std::string size3() const { return "N"; }

            void kernel_arguments(std::string & arguments_string) const{
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

        void init_rhs_global_ptr(utils::kernel_generation_stream &  stream, detail::mapped_matrix const & mat, unsigned int ks_rhs,unsigned int ns_rhs, unsigned int nl_rhs) {
          if(mat.is_row_major())
            for(unsigned int k = 0 ; k < ks_rhs ; ++k){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(k);
              stream << "__global " << mat.scalartype() << profile_.vectorization_ << " * " << ptr_name << " = " << mat.name() << " + " ;
              if(mat.is_transposed()){
                std::string i = utils::to_string(k) + " + " + offset_n + " +  get_group_id(1)*" + utils::to_string(nl_rhs);
                stream << mat.offset(i, "0", profile::size2(), profile::size3());
              }
              else{
                std::string i = utils::to_string(k);
                std::string j = offset_n + " +  get_group_id(1)*" + utils::to_string(nl_rhs);
                stream << mat.offset(i, j, profile::size2(), profile::size3());
              }
              stream << ";" << std::endl;
              //            mat->private_value(k,"*" + ptr_name);
            }
          else
            for(unsigned int n = 0 ; n < ns_rhs ; ++n){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(n);
              stream << "__global " << mat.scalartype() << profile_.vectorization_ << " * " << ptr_name << " = " << mat.name() << " +  " ;
              if(mat.is_transposed()){
                std::string i = offset_n + " +  get_group_id(1)*" + utils::to_string(nl_rhs);
                std::string j = utils::to_string(n);
                stream << mat.offset(i, j,profile::size2(), profile::size3());
              }
              else{
                std::string j = offset_n + " +  get_group_id(1)*" + utils::to_string(nl_rhs) + " + " + utils::to_string(n);
                stream << mat.offset("0",j,profile::size2(), profile::size3());
              }
              stream << ";" << std::endl;
              //            mat->private_value(n,"*" + ptr_name);
            }
        }

        void update_rhs_global_ptr(utils::kernel_generation_stream & stream, detail::mapped_matrix const & mat, unsigned int ks, unsigned int ns_rhs, unsigned int ks_rhs){
          if(mat.is_row_major() && !mat.is_transposed())
            for(unsigned int k=0 ; k<ks ; ++k)
              stream << mat.name() << "_ptr_" << k << " += " << ks_rhs << "*" << profile::size3() << " - " << ns_rhs << ";" << std::endl;
          else if(mat.is_transposed() && !mat.is_row_major())
            for(unsigned int n=0 ; n<ns_rhs ; ++n)
              stream << mat.name() << "_ptr_" << n << " += " << ns_rhs << "*" << profile::size2() << " - " << ks_rhs << ";" << std::endl;
        }


        void init_lhs_global_ptr(utils::kernel_generation_stream & stream, detail::mapped_matrix const & mat, unsigned int ms_lhs, unsigned int ks_lhs, unsigned int ml_lhs) {
          if(mat.is_row_major()){
            for(unsigned int m=0; m<ms_lhs; ++m){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(m);
              stream << "__global " << aligned_scalartype << " * " << ptr_name << " = " << mat.name() << " + ";
              if(mat.is_transposed()){
                std::string i = utils::to_string(m);
                std::string j = "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset_m;
                stream << mat.offset(i,j,profile::size1(),profile::size2());
              }
              else{
                std::string i = "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset_m + "+" + utils::to_string(m);
                stream << mat.offset(i,"0",profile::size1(),profile::size2());
              }
              stream << ";" << std::endl;
              //            mat.private_value(m,"*" + ptr_name);
            }
          }
          else{
            for(unsigned int k=0; k<ks_lhs; ++k){
              std::string ptr_name = mat.name() + "_ptr_" + utils::to_string(k);
              stream << "__global " << aligned_scalartype << " * " << ptr_name << " = " << mat.name() << " + " ;
              if(mat.is_transposed()){
                std::string j = utils::to_string(k) + "+" + "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset_m ;
                stream << mat.offset("0",j,profile::size1(),profile::size2());
              }
              else{
                std::string i = "get_group_id(0)*" + utils::to_string(ml_lhs) + "+" + offset_m;
                std::string j = utils::to_string(k);
                stream << mat.offset(i,j,profile::size1(),profile::size2());
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
                                std::string const & lmem_size2,
                                std::string const & offset,
                                unsigned int bound1,
                                unsigned int bound2,
                                symbolic_expression_tree_base const & mat_expression,
                                MatContainerT & matrices){
          std::string aligned_scalartype = (*matrices.begin())->aligned_scalartype();
          std::string scalartype = (*matrices.begin())->scalartype();
          std::string internal_size2 = (*matrices.begin())->internal_size2();
          std::string internal_size1 = (*matrices.begin())->internal_size1();
          kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          kss << "for(unsigned int i = get_local_id(0)" << " ; i < " << bound1 << "; i+= get_local_size(0)){" << std::endl;
          kss.inc_tab();
          kss << "for(unsigned int j = get_local_id(1)" << " ; j < " << bound2 << "; j+= get_local_size(1)){" << std::endl;
          kss.inc_tab();
          if((*matrices.begin())->is_rowmajor()){
            for(typename MatContainerT::iterator it = matrices.begin() ; it!=matrices.end(); ++it){
              (*it)->private_value(0,(*it)->name() +  "[" + offset + " + j  + " + internal_size2 + "*i]");
            }
            kss << aligned_scalartype << " val = " << mat_expression.generate(0) << ";" << std::endl;
            kss << "__local " << scalartype << "* ptr = " << lmem_name << " + i*" << lmem_size2 << "+j*" << profile_.vectorization_<<";" <<std::endl;
            for(unsigned int a = 0 ; a < profile_.vectorization_ ; ++a){
              if(profile_.vectorization_>1)
                kss << "*ptr++ =  val.s" << a << ";" << std::endl;
              else
                kss << "*ptr++ =  val;" << std::endl;
            }
          }
          else{
            for(typename MatContainerT::iterator it = matrices.begin() ; it!=matrices.end(); ++it){
              (*it)->private_value(0,(*it)->name() + "[" + offset + "+ j*" + internal_size1 + " + i]");
            }
            kss << aligned_scalartype << " val = " << mat_expression.generate(0) << ";" << std::endl;
            kss << "__local " << scalartype << "* ptr = " << lmem_name << " + i*" << profile_.vectorization_ * lmem_size2 << "+ j;" <<std::endl;
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
        matrix_product( const & s, profile const & p) : template_base(s, profile_), profile_(p){ }

        void core(utils::kernel_generation_stream& stream) const{

        }

      private:
        profile profile_;
    };

  }

#endif
