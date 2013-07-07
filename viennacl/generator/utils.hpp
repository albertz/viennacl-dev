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
              kgenstream(std::ostream& final_destination
                         ,unsigned int const & tab_count) : final_destination_(final_destination)
              ,tab_count_(tab_count){ }
              ~kgenstream() {  pubsync(); }
              int sync() {
                for(unsigned int i=0 ; i<tab_count_;++i)
                  final_destination_ << '\t';
                final_destination_ << str();
                str("");
                return !final_destination_;
              }
            private:
              std::ostream& final_destination_;
              unsigned int const & tab_count_;
          };

        public:
          kernel_generation_stream(std::ostream& final_destination) : std::ostream(new kgenstream(final_destination,tab_count_))
          , tab_count_(0){ }
          ~kernel_generation_stream(){ delete rdbuf(); }
          std::string str(){
            return static_cast<std::stringbuf*>(rdbuf())->str();
          }

          void inc_tab(){ ++tab_count_; }
          void dec_tab(){ --tab_count_; }

        private:
          unsigned int tab_count_;
      };

    }

  }

}
#endif
