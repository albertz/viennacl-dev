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


/** @file viennacl/generator/matrix_product.hpp
 *
 * Kernel template for the matrix product operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/mapped_objects.hpp"
#include "viennacl/generator/utils.hpp"
#include "viennacl/generator/fetch.hpp"
#include "viennacl/generator/expression_generation.hpp"
#include "viennacl/forwards.h"

#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    class matrix_product : public profile_base{

        enum access_flow{
          REGULAR,
          STRIDED
        };

        bool is_slow_impl(viennacl::ocl::device const &) const { return false; }

        std::size_t lmem_used(std::size_t scalartype_size) const {
          std::size_t lmem_used = 0;
          if(use_lhs_shared_)
            lmem_used += (ml_ + 1) * (cache_width_ + 1) * scalartype_size;
          if(use_rhs_shared_)
            lmem_used += (cache_width_ + 1) * (nl_ + 1) * scalartype_size;
          return lmem_used;
        }

        virtual void print(std::ostream & s) const{
          s << "{vector_type, local_size1, cache_width, local_size2, ms, ks, ns, use_lhs_shared, use_rhs_shared} = {"
            << simd_width_ << ","
            << local_size1_ << ", "
            << cache_width_ << ", "
            << local_size2_ << ", "
            << ms_ << ", "
            << ks_ << ", "
            << ns_ << ", "
            << use_lhs_shared_ << ", " << use_rhs_shared_ << "}" ;
        }


        bool invalid_impl(viennacl::ocl::device const & /*dev*/, size_t /*scalartype_size*/) const{
          static const unsigned int alignment = 128;
          return ml_ > alignment
              || cache_width_ > alignment
              || nl_ > alignment
              || ml_ < ms_
              || cache_width_ < ks_
              || nl_ < ns_
              || (ms_ % simd_width_) > 0
              || (ks_ % simd_width_) > 0
              || (ns_ % simd_width_) > 0;
        }

      public:
        /** @brief The user constructor */
        matrix_product(unsigned int vectorization
                , std::size_t local_size1, std::size_t cache_width, std::size_t local_size2
                , unsigned int ms, unsigned int ks, unsigned int ns
                , bool use_lhs_shared, bool use_rhs_shared) : profile_base(vectorization,local_size1, local_size2,1){
          local_size1_ = local_size1;
          local_size2_ = local_size2;
          cache_width_=cache_width;
          ml_= ms*local_size1;
          nl_=ns*local_size2;
          ms_ = ms;
          ks_=ks;
          ns_=ns;
          use_lhs_shared_ = use_lhs_shared;
          use_rhs_shared_ = use_rhs_shared;
        }

        static std::string csv_format() {
          return "Vec,LSize1,CacheWidth,LSize2,mS,kS,nS,NumGroups";
        }

        std::string csv_representation() const{
          std::ostringstream oss;
          oss << simd_width_
              << "," << local_size1_
              << "," << cache_width_
              << "," << local_size2_
              << "," << ms_
              << "," << ks_
              << "," << ns_
              << "," << use_lhs_shared_
              << "," << use_rhs_shared_;
          return oss.str();
        }

        void configure_range_enqueue_arguments(std::size_t kernel_id, statements_type  const & statements, viennacl::ocl::kernel & k, unsigned int & n_arg)  const {
          //set M, N
          scheduler::statement_node const & first_node = statements.front().second;
          unsigned int M = utils::call_on_matrix(first_node.lhs, utils::internal_size1_fun());
          unsigned int N = utils::call_on_matrix(first_node.lhs, utils::internal_size2_fun());

          //set ND range
          configure_local_sizes(k, kernel_id);
          k.global_work_size(0, M/ms_);
          k.global_work_size(1, N/ns_);

          //set arguments
          //M,N
          k.arg(n_arg++, cl_uint(M));
          k.arg(n_arg++, cl_uint(N));

          //K
          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            scheduler::statement::container_type exprs = it->first.array();
            for(scheduler::statement::container_type::iterator iit = exprs.begin() ; iit != exprs.end() ; ++iit){
              if(iit->op.type==scheduler::OPERATION_BINARY_MAT_MAT_PROD_TYPE){
                scheduler::statement_node const * current_node = &(*iit);
                //The LHS of the prod is a matrix
                if(current_node->lhs.type_family==scheduler::MATRIX_TYPE_FAMILY)
                {
                  k.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size2_fun())));
                }
                else{
                  //The LHS of the prod is a matrix expression
                  current_node = &exprs[current_node->lhs.node_index];
                  if(current_node->lhs.type_family==scheduler::MATRIX_TYPE_FAMILY)
                  {
                    if(current_node->op.type==scheduler::OPERATION_UNARY_TRANS_TYPE)
                      k.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size1_fun())));
                    else
                      k.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size2_fun())));
                  }
                  else{
                    assert(false && bool("unexpected expression tree"));
                  }
                }
                return;
              }
            }
          }

        }

        static std::string size1() { return "M";  }
        static std::string size2() { return "K"; }
        static std::string size3() { return "N"; }

        void add_kernel_arguments(statements_type  const & /*statements*/, std::string & arguments_string) const{
          arguments_string += generate_value_kernel_argument("unsigned int", "M");
          arguments_string += generate_value_kernel_argument("unsigned int", "N");
          arguments_string += generate_value_kernel_argument("unsigned int", "K");
        }

      private:

        void transform_block(detail::mapped_matrix const & /*mat_infos*/, bool store_shared
                             , unsigned int & large_block_1, unsigned int & large_block_2
                             , unsigned int & small_block_1, unsigned int & small_block_2
                             , access_flow flow) const {
          if(flow==REGULAR){
            large_block_2/=simd_width_;
            if(!store_shared)
              small_block_2/=simd_width_;
          }
          else{
            large_block_1/=simd_width_;
            if(!store_shared)
              small_block_1/=simd_width_;
          }
        }


        std::string helper_variable(utils::kernel_generation_stream & stream
                                    , bool store_in_register
                                    , std::string const & type
                                    , std::string const & name
                                    , std::string const & expr) const {
          if(!store_in_register)
            return expr;
          stream << type << " " << name << " = " << expr << ";" << std::endl;
          return name;
        }

        void fetch_element_to_local_mem(utils::kernel_generation_stream & stream,
                                std::string const & lmem_name,
                                std::size_t lmem_size2,
                                std::string const & global_ptr,
                                detail::mapped_matrix const & mat,
                                access_flow flow,
                                std::string const & i,
                                std::string const & j) const {

            if(flow==REGULAR){
                stream << "val = *(" << global_ptr << " + " << j << " + " << mat.size2()  << "*" << i << ");" << std::endl;
              for(unsigned int a = 0 ; a < simd_width_ ; ++a)
                  if(simd_width_>1)
                      stream << lmem_name << "[" << i << "*" << lmem_size2 << " + " << j << "*" << simd_width_<<" + " << a << "] = val.s" << a << ";" <<std::endl;
                  else
                      stream << lmem_name << "[" << i << "*" << lmem_size2 << " + " << j << "*" << simd_width_ << "] = val" << ";" <<std::endl;
            }
            else{
              stream << "val = *(" << global_ptr << "+ " << j << "*" << mat.size1() << " + " << i << ");" << std::endl;
              for(unsigned int a = 0 ; a < simd_width_ ; ++a)
                  if(simd_width_>1)
                      stream << lmem_name << "[" << i << "*" << simd_width_*lmem_size2 << " + " << j << " + " << a*lmem_size2 << "] = val.s" << a << ";" <<std::endl;
                  else
                      stream << lmem_name << "[" << i << "*" << simd_width_*lmem_size2 << " + " << j << "] = val" << ";" <<std::endl;
            }
        }
        void fetch_to_local_mem(utils::kernel_generation_stream & stream,
                                std::string const & lmem_name,
                                std::size_t lmem_size2,
                                std::string const & global_ptr,
                                unsigned int bound1,
                                unsigned int bound2,
                                detail::mapped_matrix const & mat,
                                access_flow flow) const {
          stream << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
          stream << "{" << std::endl;
          stream << mat.simd_scalartype() << " val;" << std::endl;
          //Can unroll
          if(bound2%local_size2_==0 && bound1%local_size1_==0){
              for(unsigned int j = 0 ; j < bound2 ; j+=local_size2_){
                  for(unsigned int i = 0 ; i < bound1 ; i+=local_size1_){
                      std::string indi = "(get_local_id(0) + " + utils::to_string(i)+")";
                      std::string indj = "(get_local_id(1) + " + utils::to_string(j)+")";
                      fetch_element_to_local_mem(stream,lmem_name,lmem_size2,global_ptr,mat,flow,indi,indj);
                  }
              }
          }
          else{
              stream << "for(unsigned int j = get_local_id(1)" << " ; j < " << bound2 << "; j+= " << local_size2_ << "){" << std::endl;
              stream.inc_tab();
              stream << "for(unsigned int i = get_local_id(0)" << " ; i < " << bound1 << "; i+= " << local_size1_ << "){" << std::endl;
              stream.inc_tab();
              fetch_element_to_local_mem(stream,lmem_name,lmem_size2,global_ptr,mat,flow,"i","j");
              stream.dec_tab();
              stream << "}" << std::endl;
              stream.dec_tab();
              stream << "}" << std::endl;

          }
          stream << "}" << std::endl;
          stream << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;

        }

        void core(std::size_t /*kernel_id*/, utils::kernel_generation_stream& stream, statements_type const & statements, std::vector<detail::mapping_type> const & mapping) const {

          //////////////////
          /// INIT
          /// //////////////

          detail::mapped_matrix const * assigned = static_cast<detail::mapped_matrix const *>(mapping.at(0).at(std::make_pair(&statements.front().second,detail::LHS_NODE_TYPE)).get());
          detail::mapped_matrix_product* prod = NULL;
          detail::mapped_matrix const * lhs = NULL;
          detail::mapped_matrix const * rhs = NULL;

          bool is_lhs_transposed = false;
          bool is_rhs_transposed = false;

          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            scheduler::statement::container_type const & exprs = it->first.array();
            std::size_t i = std::distance(statements.begin(), it);
            for(scheduler::statement::container_type::const_iterator iit = exprs.begin() ; iit != exprs.end() ; ++iit){
              if(iit->op.type==scheduler::OPERATION_BINARY_MAT_MAT_PROD_TYPE){
                prod = (detail::mapped_matrix_product *)mapping.at(i).at(std::make_pair(&(*iit), detail::PARENT_NODE_TYPE)).get();
                if(iit->lhs.type_family == scheduler::COMPOSITE_OPERATION_FAMILY){
                  is_lhs_transposed = true;
                  lhs = (detail::mapped_matrix const *)mapping.at(i).at(std::make_pair(&exprs[iit->lhs.node_index],detail::LHS_NODE_TYPE)).get();
                }
                else{
                  is_lhs_transposed = false;
                  lhs = (detail::mapped_matrix const *)mapping.at(i).at(std::make_pair(&(*iit), detail::LHS_NODE_TYPE)).get();
                }

                if(iit->rhs.type_family == scheduler::COMPOSITE_OPERATION_FAMILY){
                  is_rhs_transposed = true;
                  rhs = (detail::mapped_matrix const *)mapping.at(i).at(std::make_pair(&exprs[iit->rhs.node_index], detail::LHS_NODE_TYPE)).get();
                }
                else{
                  is_rhs_transposed = false;
                  rhs = (detail::mapped_matrix const *)mapping.at(i).at(std::make_pair(&(*iit),detail::RHS_NODE_TYPE)).get();
                }

              }
            }
          }

          if(simd_width_>1){
            std::string StrV = "/"+utils::to_string(simd_width_) ;

            for(detail::mapping_type::const_iterator it = mapping.front().begin() ; it != mapping.front().end() ; ++it){
              if(detail::mapped_matrix const * p = dynamic_cast<detail::mapped_matrix const *>(it->second.get())){
                if(p->is_row_major())
                  p->bind_sizes("M", "N"+StrV);
                else
                  p->bind_sizes("M"+StrV, "N");
              }
            }

            if(lhs->is_row_major())
              if(is_lhs_transposed)
                lhs->bind_sizes("M"+StrV, "K");
              else
                lhs->bind_sizes("M", "K"+StrV);
            else
              if(is_lhs_transposed)
                lhs->bind_sizes("M", "K"+StrV);
              else
                lhs->bind_sizes("M"+StrV, "K");


            if(rhs->is_row_major())
              if(is_rhs_transposed)
                rhs->bind_sizes("K"+StrV, "N");
              else
                rhs->bind_sizes("K", "N"+StrV);
            else
              if(is_rhs_transposed)
                rhs->bind_sizes("K", "N"+StrV);
              else
                rhs->bind_sizes("K"+StrV, "N");


          }
          else{
            for(detail::mapping_type::const_iterator it = mapping.front().begin() ; it != mapping.front().end() ; ++it){
              if(detail::mapped_matrix const * p = dynamic_cast<detail::mapped_matrix const *>(it->second.get())){
                p->bind_sizes("M", "N");
              }
            }

            lhs->bind_sizes("M", "K");
            rhs->bind_sizes("K", "N");
          }




          access_flow result_access_flow;
          if(assigned->is_row_major())
            result_access_flow = REGULAR;
          else
            result_access_flow = STRIDED;

          access_flow lhs_access_flow;
          if((lhs->is_row_major() && !is_lhs_transposed)
             ||(!lhs->is_row_major() && is_lhs_transposed))
            lhs_access_flow = REGULAR;
          else
            lhs_access_flow = STRIDED;

          access_flow rhs_access_flow;
          if((rhs->is_row_major() && !is_rhs_transposed)
             ||(!rhs->is_row_major() && is_rhs_transposed))
            rhs_access_flow = REGULAR;
          else
            rhs_access_flow = STRIDED;


          unsigned int ml_res = ml_, nl_res = nl_, ms_res = ms_, ns_res = ns_;
          unsigned int ml_lhs = ml_, cache_width_lhs = cache_width_, ms_lhs = ms_, ks_lhs = ks_;
          unsigned int cache_width_rhs = cache_width_, nl_rhs = nl_, ks_rhs = ks_, ns_rhs = ns_;

          transform_block(*assigned,false,ml_res,nl_res,ms_res,ns_res,result_access_flow);
          transform_block(*lhs,use_lhs_shared_,ml_lhs,cache_width_lhs,ms_lhs,ks_lhs,lhs_access_flow);
          transform_block(*rhs,use_rhs_shared_,cache_width_rhs,nl_rhs,ks_rhs,ns_rhs,rhs_access_flow);

          //////////////////
          /// DECLARATIONS
          /// //////////////


          std::size_t local_lhs_size1 = ml_ ;
          std::size_t local_lhs_size2 = cache_width_ + 1;

          std::size_t local_rhs_size1 = cache_width_;
          std::size_t local_rhs_size2 = nl_ + 1;

          ///Result Values
          for(unsigned int m=0; m< ms_res; ++m)
            for(unsigned int n=0; n < ns_res ; ++n)
              stream << assigned->simd_scalartype() << " " << "res" << m << "_" << n << " = (" << assigned->simd_scalartype() << ")(0) ;" << std::endl;

          ///Local memory
          if(use_lhs_shared_)
            stream << "__local " << lhs->scalartype() << " lhs_buf[" << local_lhs_size1*local_lhs_size2 << "]" << ";" << std::endl;
          if(use_rhs_shared_)
            stream << "__local " << rhs->scalartype() << " rhs_buf[" << local_rhs_size1*local_rhs_size2 << "]" << ";" << std::endl;

          ///LHS - Local Memory Offset
          if(use_lhs_shared_){
            std::string i = "get_group_id(0)*" + utils::to_string(ml_lhs);
            stream << "__global " << lhs->simd_scalartype() << "* global_lhs_ptr = " << lhs->name() << " + ";
            if(lhs_access_flow==REGULAR)
              stream << "(" << i << ")" << "*" << lhs->size2();
            else
              stream << i;
            stream << ";" << std::endl;
          }

          ///LHS - Global Memory pointer
          else{
            if(lhs_access_flow==REGULAR)
              for(unsigned int m=0; m<ms_lhs; ++m)
                stream << "__global " << lhs->simd_scalartype() << "* " << "lhs_ptr_" << m << " = " << lhs->name() << " + "
                       << lhs->size2() << "* ("
                       << "get_group_id(0)*" << ml_lhs << "+" << "get_local_id(0)*" << ms_lhs << "+" << m
                       << " );" << std::endl;
            else
              for(unsigned int k=0; k<ks_lhs; ++k)
                stream << "__global " << lhs->simd_scalartype() << "* " << "lhs_ptr_" << k << " = " << lhs->name() << " + "
                       << "(" << lhs->size1() << ")*" << k
                       << "+ " << "get_group_id(0)*" << ml_lhs << "+" << "get_local_id(0)*" << ms_lhs << ";" << std::endl;
          }

          ///RHS - Local Memory Offset
          if(use_rhs_shared_){
            std::string j = "get_group_id(1)*" + utils::to_string(nl_rhs);
            stream << "__global " << rhs->simd_scalartype() << "* global_rhs_ptr = " << rhs->name() << " + ";
            if(rhs_access_flow==REGULAR)
              stream << j;
            else
              stream << "(" << j << ")" << "*" << rhs->size1();
            stream << ";" << std::endl;
          }

          ///RHS - Global Memory Pointer
          else{
            if(rhs_access_flow==REGULAR)
              for(unsigned int k = 0 ; k < ks_rhs ; ++k)
                stream << "__global " << rhs->simd_scalartype() << "* " << "rhs_ptr_" << k << " = " << rhs->name() << " + "
                       << "(" << k << ")" << "*" << rhs->size2()
                       << " + " << "get_local_id(1)*" << ns_rhs << " + get_group_id(1)*" << nl_rhs
                       << ";" << std::endl;
            else
              for(unsigned int n = 0 ; n < ns_rhs ; ++n)
                stream << "__global " << rhs->simd_scalartype() << "* " << "rhs_ptr_" << n << " = " << rhs->name() << " +  "
                       << "(" << "get_local_id(1)*" << ns_rhs << " + get_group_id(1)*" << nl_rhs << " + " << n << ")" << "*" << rhs->size1()
                       << ";" << std::endl;
          }


          ///Large Work-group Wise loop
          std::string block_num = helper_variable(stream,false,"unsigned int", "block_num", "K/" + utils::to_string(cache_width_));
          stream << "for(unsigned int bl=0 ; bl<" << block_num << " ; ++bl){" << std::endl;
          stream.inc_tab();

          ///Update LHS Local Memory and pointers (if necessary)
          if(use_lhs_shared_){
            fetch_to_local_mem(stream,"lhs_buf",local_lhs_size2,"global_lhs_ptr",ml_lhs,cache_width_lhs,*lhs,lhs_access_flow);
            for(unsigned int m=0; m<ms_lhs; ++m)
              stream << "__local " << lhs->scalartype() << "* lhs_ptr_" << m << " = lhs_buf + "
                     << "(" << "get_local_id(0)*" << ms_lhs << "+" << m << ")" << "*" << local_lhs_size2
                     << ";" << std::endl;
          }

          ///Update RHS Local Memory and pointers (if necessary)
          if(use_rhs_shared_){
            fetch_to_local_mem(stream,"rhs_buf", local_rhs_size2, "global_rhs_ptr",cache_width_rhs,nl_rhs,*rhs,rhs_access_flow);
            for(unsigned int k=0; k<ks_rhs; ++k)
              stream << "__local " << rhs->scalartype() << "* rhs_ptr_" << k << " = rhs_buf + "
                     << k*local_rhs_size2 << " + " << "get_local_id(1)*" << ns_rhs
                     << ";" << std::endl;
          }


          stream << " for(unsigned int bs=0 ; bs < " << cache_width_/ks_  << " ; ++bs){" << std::endl;
          stream.inc_tab();


          for(unsigned int k = 0 ; k < ks_rhs ; ++k){
            for(unsigned int n=0 ; n < ns_rhs ; ++n){
              if(use_rhs_shared_ )
                  stream << rhs->scalartype() << " val_rhs_" << k << "_" << n << " = * rhs_ptr_" << k << "++";
              else{
                stream << rhs->simd_scalartype() << " val_rhs_" << k << "_" << n << " = " ;
                if(rhs_access_flow==REGULAR)
                  stream << "* rhs_ptr_" << k << "++";
                else
                  stream  << "* rhs_ptr_" << n << "++";
              }
              stream << ";";
              stream << std::endl;
            }
          }


          for(unsigned int k = 0 ; k < ks_lhs ; ++k){
            for(unsigned int m=0 ; m < ms_lhs ; ++m){
              if(use_lhs_shared_)
                stream << lhs->scalartype() << " " << "val_lhs_" << m << "_" << k << " = * lhs_ptr_" << m << "++" ;
              else{
                  stream << lhs->simd_scalartype() << " " << "val_lhs_" << m << "_" << k << " = ";
                  if(lhs_access_flow==REGULAR)
                    stream << "* lhs_ptr_" << m << "++";
                  else
                    stream << "* lhs_ptr_" << k << "++";
              }
              stream << ";";
              stream << std::endl;
            }
          }

          for(unsigned int k = 0 ; k < ks_ ; ++k){
              for(unsigned int m=0 ; m < ms_ ; ++m){
                  for(unsigned int n=0 ; n < ns_ ; ++n){
                      std::ostringstream res_oss;
                      if(simd_width_==1)
                          res_oss << "res" << m << "_" << n ;
                      else{
                          if(result_access_flow==REGULAR)
                              res_oss << "res" << m << "_" << n/simd_width_  << ".s" << n%simd_width_;
                          else
                              res_oss << "res" << m/simd_width_ << "_" << n  << ".s" << m%simd_width_;
                      }

                      std::ostringstream lhs_oss;
                      if(use_lhs_shared_ || simd_width_==1){
                          lhs_oss << "val_lhs_" << m << "_" << k;
                      }
                      else{
                          if(lhs_access_flow==REGULAR)
                              lhs_oss << "val_lhs_" << m << "_" << k/simd_width_ << ".s" << k%simd_width_;
                          else
                              lhs_oss << "val_lhs_" << m/simd_width_ << "_" << k << ".s" << m%simd_width_;
                      }

                      std::ostringstream rhs_oss;
                      if(use_rhs_shared_ || simd_width_==1){
                          rhs_oss << "val_rhs_" << k << "_" << n;
                      }
                      else{
                          if(rhs_access_flow==REGULAR)
                              rhs_oss << "val_rhs_" << k << "_" << n/simd_width_ << ".s" << n%simd_width_;
                          else
                              rhs_oss << "val_rhs_" << k/simd_width_ << "_" << n << ".s" << k%simd_width_;
                      }


                      stream << res_oss.str() << "+=" << lhs_oss.str() << "*" << rhs_oss.str() << ";" << std::endl;
                  }
              }
          }


          if(use_rhs_shared_){
            for(unsigned int k=0 ; k<ks_ ; ++k)
              stream << "rhs_ptr_" << k << " += " << ks_rhs*local_rhs_size2 - ns_rhs << ";" << std::endl;
          }
          else{
            if(rhs_access_flow==REGULAR)
              for(unsigned int k=0 ; k<ks_ ; ++k)
                stream << "rhs_ptr_" << k << " += " << ks_rhs << "*" << rhs->size2() << " - " << ns_rhs << ";" << std::endl;
          }

          if(!use_lhs_shared_){
            if(lhs_access_flow==STRIDED)
              for(unsigned int k=0 ; k<ks_lhs ; ++k)
                stream << "lhs_ptr_" << k << " += " << ks_lhs << "*" << lhs->size1() << " - " << ms_lhs << ";" << std::endl;
          }



          stream.dec_tab();
          stream << "}" << std::endl;

          if(use_lhs_shared_){
            if(lhs_access_flow==REGULAR)
              stream << "global_lhs_ptr += " << cache_width_lhs << ";" << std::endl;
            else
              stream << "global_lhs_ptr += " << cache_width_lhs << "*" << lhs->size1() << ";" << std::endl;
          }

          if(use_rhs_shared_){
            if(rhs_access_flow==REGULAR)
              stream << "global_rhs_ptr += " << cache_width_rhs << "*" << rhs->size2() << ";" << std::endl;
            else
              stream << "global_rhs_ptr += " << cache_width_rhs << ";" << std::endl;
          }

          stream.dec_tab();
          stream << "}" << std::endl;

          for(unsigned int m=0 ; m < ms_res ; ++m){
            for(unsigned int n=0 ; n < ns_res ; ++n){
              std::string i = "get_global_id(0)*" + utils::to_string(ms_res) + "+" + utils::to_string(m);
              std::string j = "get_global_id(1)*" + utils::to_string(ns_res) + "+" + utils::to_string(n);
              prod->access_name("res"+utils::to_string(m)+"_"+utils::to_string(n));
              std::string str;
              detail::traverse(statements.front().first, statements.front().second, detail::expression_generation_traversal(std::make_pair(i, j), -1, str, mapping[0]), false);
              stream << str << ";" << std::endl;
            }
          }


        }

      private:
        std::size_t local_size1_;
        std::size_t local_size2_;
        std::size_t cache_width_;

        std::size_t ml_;
        std::size_t nl_;

        std::size_t ms_;
        std::size_t ks_;
        std::size_t ns_;

        bool use_lhs_shared_;
        bool use_rhs_shared_;
    };

  }

}

#endif
