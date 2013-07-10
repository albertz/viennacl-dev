#ifndef VIENNACL_META_IS_BASE_OF_HPP_
#define VIENNACL_META_IS_BASE_OF_HPP_

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

/** @file is_base_of.hpp
    @brief Implementation of is_base_of for compile-time relationship checking
*/

namespace viennacl
{
  template <typename B, typename D>
  struct is_base_of
  {
    private:
      typedef char (&yes)[1];
      typedef char (&no)[2];

      template <typename B_, typename D_>
      struct Host
      {
        operator B_*() const;
        operator D_*();
      };

    public:
      template <typename T>
      static yes check(D*, T);
      static no check(B*, int);

     static const bool value = sizeof(check(Host<B,D>(), int())) == sizeof(yes);
  };

} //namespace viennacl


#endif
