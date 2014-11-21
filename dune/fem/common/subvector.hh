#ifndef DUNE_FEM_COMMON_SUBVECTOR_HH
#define DUNE_FEM_COMMON_SUBVECTOR_HH

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/densevector.hh>
#include <dune/common/nullptr.hh>

namespace Dune
{

  namespace Fem
  {
    template< class V >
    class DenseSubVector;
  }


  template< class V >
  struct DenseMatVecTraits< Fem::DenseSubVector< V > >
  {
    typedef Fem::DenseSubVector< V > derived_type;
    typedef V container_type;

    typedef typename V::value_type value_type;
    typedef typename V::size_type size_type;
  };

  template< class V >
  struct FieldTraits< Fem::DenseSubVector< V > >
  {
    typedef typename FieldTraits< typename V::value_type >::field_type field_type;
    typedef typename FieldTraits< typename V::value_type >::real_type real_type;
  };

  namespace Fem
  {

    // DenseSubVector
    // ---------

    template< class V >
    class DenseSubVector : public DenseVector< DenseSubVector< V > >
    {
      typedef DenseSubVector< V > This;
      typedef DenseVector< DenseSubVector < V > > Base;


    public:
      typedef typename Base::size_type size_type;
      typedef typename Base::value_type value_type;

      //! Constructor making uninitialized vector
      explicit DenseSubVector ( V& v, size_type size, size_type offset )
      : v_( v ),
        size_( size ),
        offset_( offset )
      {}

      DenseSubVector ( const This &other )
      : v_( other.v_ ),
        size_( other.size_ ),
        offset_( other.offset_ )
      {}

      using Base::operator=;

      //==== forward some methods of std::vector
      /** \brief Number of elements for which memory has been allocated.

          capacity() is always greater than or equal to size().
       */

      void resize ( size_type ) {}

      //==== make this thing a vector
      size_type vec_size () const { return size_; }

      value_type &vec_access ( size_type i ) { return v_[ i + offset_ ]; }
      const value_type &vec_access ( size_type i ) const { return v_[ i + offset_ ]; }
    private:
      typename DenseMatVecTraits< This > :: container_type &v_;
      const size_type size_, offset_;
    };

    /** @} end documentation */

  } // namespace Fem


} // namespace Dune

#endif //#ifndef DUNE_FEM_COMMON_SUBVECTOR_HH
