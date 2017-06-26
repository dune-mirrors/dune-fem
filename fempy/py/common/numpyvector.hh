#ifndef DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH
#define DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fempy/pybind11/pybind11.hh>

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
        : array_( pybind11::buffer_info( nullptr, sizeof( T ),
                  pybind11::format_descriptor< T >::value, 1, { size }, { sizeof( T ) } )
                ),
          // bufferInfo_( array_.request( true ) ),
          size_(size)
      {}

      NumPyVector ( pybind11::buffer buf )
        : array_( buf ), // buf.request(true) ),
          size_( 0 )
      {
        pybind11::buffer_info info = buf.request();
        if (info.ndim != 1)
          DUNE_THROW( InvalidStateException, "NumPyVector can only be created from one-dimensional array" );
        size_ = info.shape[0];
      }

      NumPyVector ( const This &other ) = delete;
      NumPyVector ( This &&other ) = delete;

      ~NumPyVector()
      {
        std::cout << "in numpy destructor" << std::endl;
      }

      This &operator= ( const This &other ) = delete;
      This &operator= ( This &&other ) = delete;

      operator pybind11::array_t< T > () const { return array_; }

      const value_type &operator[] ( size_type index ) const
      {
        return data()[ index ];
      }
      value_type &operator[] ( size_type index )
      {
        return data()[ index ];
      }
      value_type &vec_access ( size_type index )
      {
        return data()[ index ];
      }
      const value_type &vec_access ( size_type index ) const
      {
        return data()[ index ];
      }

      inline const value_type *data () const
      {
        return static_cast< value_type * >( array_.request(false).ptr ); // bufferInfo_.ptr );
      }
      inline value_type *data ()
      {
        return static_cast< value_type * >( array_.request(true).ptr ); // bufferInfo_.ptr );
      }
      pybind11::array_t< T > &coefficients()
      {
        return array_;
      }
      pybind11::array_t< T > &coefficients() const
      {
        return array_;
      }

      void reserve ( size_type newCapacity ) {}

      void resize ( size_type newSize, value_type value = value_type() )
      {
        DUNE_THROW(NotImplemented, "resize is not available on a NumpyVector at the moment");
        if( newSize == size() )
          return;

        pybind11::array_t< T > other( pybind11::buffer_info( nullptr, sizeof( T ),
                      pybind11::format_descriptor< T >::value, 1, { newSize }, { sizeof( T ) } ) );
        value_type* buffer = static_cast< value_type * >( array_.request(true).ptr );
        for( std::size_t i = 0; i < std::min( newSize, size() ); ++i )
          buffer[ i ] = (*this)[ i ];
        for( std::size_t i = std::min( newSize, size() ); i < newSize; ++i )
          buffer[ i ] = value;
        // bufferInfo_ = other.request( true ); // causes corrupted double list error at end of program
        array_ = std::move( other );
        size_ = newSize;
        std::cout << "Sizes: " << size() << " " << newSize << std::endl;
      }

      size_type size () const
      {
        return size_; // bufferInfo_.shape[ 0 ];
      }
      size_type vec_size () const
      {
        return size_; // bufferInfo_.shape[ 0 ];
      }

    private:
      mutable pybind11::array_t< T > array_;
      size_type size_;
      // pybind11::buffer_info bufferInfo_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH
