#ifndef DUNE_FEM_COMMON_REFERENCEVECTOR_HH
#define DUNE_FEM_COMMON_REFERENCEVECTOR_HH

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/densevector.hh>

namespace Dune
{

  namespace Fem
  {
    template< class K, class A = std::allocator< K* > >
    class DynamicReferenceVector;
  }


  template< class K, class A >
  struct DenseMatVecTraits< Fem::DynamicReferenceVector< K, A > >
  {
    typedef Fem::DynamicReferenceVector< K, A > derived_type;
    typedef std::vector< K*, typename A::template rebind< K* >::other > container_type;

    typedef K value_type;
    typedef typename container_type::size_type size_type;
  };

  template< class K, class A >
  struct FieldTraits< Fem::DynamicReferenceVector< K, A > >
  {
    typedef typename FieldTraits< K >::field_type field_type;
    typedef typename FieldTraits< K >::real_type real_type;
  };

  namespace Fem
  {

    /** @addtogroup DenseMatVec
        @{
     */

    /*! \file
     * \brief This file implements a dense vector with a dynamic size.
     */


    /** \brief Construct a vector with a dynamic size.
     *
     * \tparam K is the field type (use float, double, complex, etc)
     */
    template< class K, class A >
    class DynamicReferenceVector : public DenseVector< DynamicReferenceVector< K, A > >
    {
      typedef DynamicReferenceVector< K, A > This;
      typedef DenseVector< DynamicReferenceVector < K, A > > Base;

      typename DenseMatVecTraits< This > :: container_type data_;
    public:
      typedef typename Base::size_type size_type;
      typedef typename Base::value_type value_type;

      //! Constructor making uninitialized vector
      explicit DynamicReferenceVector ( const A &a = A() )
      : data_( a )
      {}

      //! Constructor making vector with identical coordinates
      explicit DynamicReferenceVector ( size_type n, const A &a = A() )
      : data_( n, nullptr, a )
      {}

      DynamicReferenceVector ( const This &other )
      : data_( other.data_ )
      {}

      DynamicReferenceVector ( This &&other )
      : data_( std::move( other.data_ ) )
      {}

      //! assignment
      using Base::operator=;

      template< class V >
      This & operator= ( const DenseVector< V > &other )
      {
        assert( data_.size() == other.size() );
        std::copy( other.begin(), other.end(), Base::begin() );
        return *this;
      }

      This & operator= ( const This &other )
      {
        assert( data_.size() == other.size() );
        std::copy( other.begin(), other.end(), Base::begin() );
        return *this;
      }

      This & operator= ( This &&other )
      {
        data_ = std::move( other.data_ );
        return *this;
      }

      //==== forward some methods of std::vector
      /** \brief Number of elements for which memory has been allocated.

          capacity() is always greater than or equal to size().
       */
      size_type capacity () const { return data_.capacity(); }

      void resize ( size_type n ) { data_.resize( n, nullptr ); }
      void reserve ( size_type n ) { data_.reserve( n ); }

      //! bind i-th entry to a reference
      void bind ( size_type i, K& u ) { assert( i < data_.size() ); data_[ i ] = &u; }
      void unbind ( size_type i ) { asssert( i < data_.size() ); data_[ i ] = nullptr; }

      //==== make this thing a vector
      size_type size () const { return data_.size(); }

      K &operator[] ( size_type i ) { return *data_[ i ]; }
      const K &operator[] ( size_type i ) const { return *data_[ i ]; }
    };

    /** @} end documentation */

  } // namespace Fem


} // namespace Dune

#endif //#ifndef DUNE_FEM_COMMON_REFERENCEVECTOR_HH
