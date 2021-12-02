#warning "This header should not be needed anymore. Remove it from the include list!"

#ifndef DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH
#define DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fempy/pybind11/pybind11.hh>
#include <dune/python/common/numpyvector.hh>

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
      : public Dune::Python::NumPyVector< T >
    {
      typedef NumPyVector< T > This;
      typedef Dune::Python::NumPyVector< T > Base;

    public:
      typedef typename Base::size_type size_type;
      typedef typename Base::value_type value_type;

      using Base::size;
      using Base::data;

      explicit NumPyVector ( size_type size ) : Base( size )
      {}

      NumPyVector ( pybind11::buffer buf ) : Base( buf )
      {}

      NumPyVector ( const This &other ) = delete;
      NumPyVector ( This &&other ) = delete;

      ~NumPyVector()
      {
        //std::cout << "in numpy destructor" << std::endl;
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
        const std::size_t minSize = std::min( newSize, size() );
        // copy old values up to min of size and newSize
        for( std::size_t i = 0; i < minSize; ++i )
          buffer[ i ] = data()[ i ];

        // initialize all other values with default
        for( std::size_t i = minSize; i < newSize; ++i )
          buffer[ i ] = value;

        // bufferInfo_ = other.request( true ); // causes corrupted double list error at end of program
        array_ = std::move( other );
        // re-capture data pointer
        dataPtr_ = static_cast< value_type * >( array_.request(true).ptr );
        // store size
        size_ = newSize;

        //std::cout << "Sizes: " << size() << " " << newSize << std::endl;
      }

    protected:
      using Base :: array_;
      using Base :: dataPtr_;
      using Base :: size_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_COMMON_NUMPYVECTOR_HH
