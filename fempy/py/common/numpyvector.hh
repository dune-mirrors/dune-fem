#ifndef DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH
#define DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fempy/pybind11/numpy.h>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class T >
    class NumPyVector;

  } // namespace FemPy



  // DenseMatVecTraits for NumPyVector
  // ---------------------------------

  template< class T >
  struct DenseMatVecTraits< FemPy::NumPyVector< T > >
  {
    typedef FemPy::NumPyVector< T > derived_type;
    typedef pybind11::array_t< T > container_type;
    typedef T value_type;
    typedef std::size_t size_type;
  };



  // FieldTraits for NumPyVector
  // ---------------------------

  template< class T >
  struct FieldTraits< FemPy::NumPyVector< T > >
  {
    typedef typename FieldTraits< T >::field_type field_type;
    typedef typename FieldTraits< T >::real_type real_type;
  };


  namespace FemPy
  {

    template< class T >
    class NumPyVector
      : public DenseVector< NumPyVector< T > >
    {
      typedef NumPyVector< T > This;
      typedef DenseVector< NumPyVector< T > > Base;

    public:
      typedef typename Base::size_type size_type;
      typedef typename Base::value_type value_type;

      explicit NumPyVector ( size_type size )
        : array_( pybind11::buffer_info( nullptr, sizeof( T ), pybind11::format_descriptor< T >::value, 1, { size }, { sizeof( T ) } ) ),
          bufferInfo_( array_.request( true ) )
      {}

      NumPyVector ( pybind11::array_t< T > array )
        : array_( std::move( array ) ),
          bufferInfo_( array_.request( true ) )
      {
        if( bufferInfo_.ndim() != 1 )
          DUNE_THROW( InvalidStateException, "NumPyVector can only be created from one-dimensional array" );
      }

      NumPyVector ( const This &other ) = delete;
      NumPyVector ( This &&other ) = default;

      This &operator= ( const This &other ) = delete;
      This &operator= ( This && ) = default;

      operator pybind11::array_t< T > () const { return array_; }

      const value_type &operator[] ( size_type index ) const
      {
        return static_cast< T * >( bufferInfo_.ptr )[ index * bufferInfo_.strides[ 0 ] ];
      }

      value_type &operator[] ( size_type index )
      {
        return static_cast< T * >( bufferInfo_.ptr )[ index * bufferInfo_.strides[ 0 ] ];
      }
      value_type &vec_access ( size_type index )
      {
        return static_cast< T * >( bufferInfo_.ptr )[ index * bufferInfo_.strides[ 0 ] ];
      }
      const value_type &vec_access ( size_type index ) const
      {
        return static_cast< T * >( bufferInfo_.ptr )[ index * bufferInfo_.strides[ 0 ] ];
      }

      void reserve ( size_type newCapacity ) {}

      void resize ( size_type newSize, value_type value = value_type() )
      {
        if( newSize == size() )
          return;

        NumPyVector other( newSize );
        for( std::size_t i = 0; i < std::min( newSize, size() ); ++i )
          other[ i ] = (*this)[ i ];
        for( std::size_t i = std::min( newSize, size() ); i < newSize; ++i )
          other[ i ] = value;
        *this = std::move( other );
      }

      size_type size () const
      {
        return bufferInfo_.shape[ 0 ];
      }
      size_type vec_size () const
      {
        return bufferInfo_.shape[ 0 ];
      }

    private:
      pybind11::array_t< T > array_;
      pybind11::buffer_info bufferInfo_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH
