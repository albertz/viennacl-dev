#ifndef VIENNACL_GENERATOR_UTILS_HPP
#define VIENNACL_GENERATOR_UTILS_HPP

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


/** @file viennacl/generator/utils.hpp
    @brief Internal utils
*/

#include <sstream>

namespace viennacl{

  namespace generator{

    namespace utils{

      template <class T>
      inline std::string to_string ( T const t )
      {
        std::stringstream ss;
        ss << t;
        return ss.str();
      }

      class kernel_generation_stream : public std::ostream{
        private:

          class kgenstream : public std::stringbuf{
            public:
              kgenstream(std::ostringstream& oss,unsigned int const & tab_count) : oss_(oss), tab_count_(tab_count){ }
              int sync() {
                for(unsigned int i=0 ; i<tab_count_;++i)
                  oss_ << '\t';
                oss_ << str();
                str("");
                return !oss_;
              }
              ~kgenstream() {  pubsync(); }
            private:
              std::ostream& oss_;
              unsigned int const & tab_count_;
          };

        public:
          kernel_generation_stream() : std::ostream(new kgenstream(oss,tab_count_)), tab_count_(0){ }

          std::string str(){ return oss.str(); }

          void inc_tab(){ ++tab_count_; }

          void dec_tab(){ --tab_count_; }

          ~kernel_generation_stream(){ delete rdbuf(); }

        private:
          unsigned int tab_count_;
          std::ostringstream oss;
      };

    }

  }

}
#endif
