#ifndef DUNE_FEM_VECTOR_HH
#define DUNE_FEM_VECTOR_HH

#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

#include <dune/common/math.hh>
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/storage/array.hh>
#include <dune/fem/io/streams/streams.hh>

/*! @addtogroup VectorClasses
    @{
*/

namespace Dune
{

  namespace Fem
  {

    template< class VT >
    struct VectorInterfaceArrayTraits
    {
      typedef typename VT::VectorType ArrayType;
      typedef typename VT::FieldType ElementType;

      typedef typename VT::ConstIteratorType ConstIteratorType;
      typedef typename VT::IteratorType IteratorType;
    };


    //! An abstract vector interface
    template< class VT >
    class VectorInterface
    : public ArrayInterface< VectorInterfaceArrayTraits< VT > >
    {
      typedef VectorInterface< VT > ThisType;
      typedef ArrayInterface< VectorInterfaceArrayTraits< VT > > BaseType;

      template< class > friend class VectorInterface;

    public:
      typedef VT Traits;

      //! type of this interface
      typedef ThisType VectorInterfaceType;

      //! type of the implementation (Barton-Nackman)
      typedef typename Traits::VectorType VectorType;

      //! field type for the vector
      typedef typename Traits::FieldType FieldType;
      typedef FieldType value_type;

      //! type of constant iterator
      typedef typename Traits::ConstIteratorType ConstIteratorType;

      //! type of iterator
      typedef typename Traits::IteratorType IteratorType;

      //! Assign another vector to this one
      template< class T >
      VectorType& operator= ( const VectorInterface< T > &v )
      {
        asImp().assign( v );
        return asImp();
      }

      //! Assign another vector to this one
      VectorType& operator= ( const ThisType &v )
      {
        asImp().assign( v );
        return asImp();
      }

      //! Initialize all fields of this vector with a scalar
      VectorType &operator= ( const FieldType &s )
      {
        asImp().assign( s );
        return asImp();
      }

      //! Returns a const reference to the field indexed by index
      const FieldType &operator[] ( unsigned int index ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
        return asImp()[ index ];
      }

      //! Returns a reference to the field indexed by index
      FieldType &operator[] ( unsigned int index )
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
        return asImp()[ index ];
      }

      //! Add another vector to this one
      template< class T >
      VectorType &operator+= ( const VectorInterface< T > &v )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator+=( v.asImp() ) );
        return asImp();
      }

      //! Subtract another vector from this one
      template< class T >
      VectorType &operator-= ( const VectorInterface< T > &v )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator-=( v.asImp() ) );
        return asImp();
      }

      //! Multiply this vector by a scalar
      VectorType &operator*= ( const FieldType &s )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator*=( s ) );
        return asImp();
      }

      //! Add a multiple of another vector to this one
      template< class T >
      VectorType &addScaled ( const FieldType &s, const VectorInterface< T > &v )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().add( s, v.asImp() ) );
        return asImp();
      }

      /** \brief copy another vector to this one
       *
       *  Copies the data from another vector to this one. Both vectors must be of
       *  the same size.
       *
       *  \param[in]  v  vector to copy
       */
      template< class T >
      void assign ( const VectorInterface< T > &v )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( v.asImp() ) );
      }

      //! Initialize all fields of this vector with a scalar
      void assign ( const FieldType &s )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( s ) );
      }

      /** \brief initialize the vector to 0 */
      void clear ()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().clear() );
      }

      //! obtain begin iterator
      ConstIteratorType begin () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      //! obtain begin iterator
      IteratorType begin ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      //! obtain end iterator
      ConstIteratorType end () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      //! obtain end iterator
      IteratorType end ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      //! Returns the vector's size
      unsigned int size () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
        return asImp().size();
      }

    protected:
      using BaseType::asImp;
    };


    template< class Vector >
    struct SupportsVectorInterface
    {
      typedef VectorInterface< typename Vector::Traits > VectorInterfaceType;
      static const bool v = std::is_convertible< Vector, VectorInterfaceType >::value;
    };


    template< class V, class W >
    struct ExtractCommonFieldType
    {
      typedef typename V::FieldType FieldType;

    private:
      static_assert( (std::is_same< FieldType, typename W::FieldType >::value),
                     "FieldType must be identical." );
    };



    template< class Field, class Vector >
    struct VectorDefaultTraits
    {
      typedef Field FieldType;
      typedef Vector VectorType;

      typedef ArrayDefaultIterator< FieldType, VectorType > IteratorType;
      typedef ArrayDefaultIterator< const FieldType, const VectorType > ConstIteratorType;
    };



    // VectorDefault
    // -------------

    /** \class VectorDefault
     *  \ingroup Vector
     *  \brief default implementation of VectorInterface
     */
    template< class Field, class Vector >
    class VectorDefault
    : public VectorInterface< VectorDefaultTraits< Field, Vector > >
    {
      typedef VectorDefault< Field, Vector > ThisType;
      typedef VectorInterface< VectorDefaultTraits< Field, Vector > > BaseType;

    public:
      typedef typename BaseType :: FieldType FieldType;

      typedef typename BaseType :: VectorInterfaceType VectorInterfaceType;
      typedef typename BaseType :: VectorType VectorType;

      typedef typename BaseType :: ConstIteratorType ConstIteratorType;
      typedef typename BaseType :: IteratorType IteratorType;

      //! Add another vector to this one
      template< class T >
      VectorType &operator+= ( const VectorInterface< T > &v )
      {
        const unsigned int size = this->size();
        assert( size == v.size() );
        for( unsigned int i = 0; i < size; ++i )
          (*this)[ i ] += v[ i ];
        return asImp();
      }

      //! Subtract another vector from this one
      template< class T >
      VectorType &operator-= ( const VectorInterface< T > &v )
      {
        const unsigned int size = this->size();
        assert( size == v.size() );
        for( unsigned int i = 0; i < size; ++i )
          (*this)[ i ] -= v[ i ];
        return asImp();
      }

      //! Multiply this vector by a scalar
      VectorType &operator*= ( const FieldType &s )
      {
        for( auto& entry : *this )
          entry *= s;
        return asImp();
      }

      //! Add a multiple of another vector to this one
      template< class T >
      VectorType &addScaled ( const FieldType &s, const VectorInterface< T > &v )
      {
        const unsigned int size = this->size();
        assert( size == v.size() );
        for( unsigned int i = 0; i < size; ++i )
          (*this)[ i ] += s * v[ i ];
        return asImp();
      }

      /** \copydoc Dune::Fem::VectorInterface::assign(const VectorInterface<T> &v) */
      template< class T >
      void assign ( const VectorInterface< T > &v )
      {
        assert( this->size() == v.size() );
        std::copy( v.begin(), v.end(), asImp().begin() );
      }

      //! Initialize all fields of this vector with a scalar
      void assign ( const FieldType &s )
      {
        std::fill( asImp().begin(), asImp().end(), s );
      }

      /** \copydoc Dune::Fem::VectorInterface::clear() */
      void clear ()
      {
        asImp().assign( 0 );
      }

      //! obtain begin iterator
      ConstIteratorType begin () const
      {
        return ConstIteratorType( asImp(), 0 );
      }

      //! obtain begin iterator
      IteratorType begin ()
      {
        return IteratorType( asImp(), 0 );
      }

      //! obtain end iterator
      ConstIteratorType end () const
      {
        return ConstIteratorType( asImp(), size() );
      }

      //! obtain end iterator
      IteratorType end ()
      {
        return IteratorType( asImp(), size() );
      }

      using BaseType :: size;

    protected:
      using BaseType :: asImp;
    };



    /** \class DynamicVector
     *  \brief A vector using a std::vector as storage
     *
     *  An implementation of VectorInterface using a std::vector to provide the
     *  fields.
     */
    template< class Field, template< class > class Allocator = std::allocator >
    class DynamicVector
    : public VectorDefault< Field, DynamicVector< Field, Allocator > >
    {
      typedef DynamicVector< Field, Allocator > ThisType;
      typedef VectorDefault< Field, ThisType > BaseType;

    public:
      //! field type of the vector
      typedef Field FieldType;

      using BaseType :: assign;

    protected:
      std::vector< FieldType, Allocator<FieldType> > fields_;

    public:
      //! Constructor setting up a vector of a specified size
      explicit DynamicVector ( unsigned int size = 0 )
      : fields_( size )
      {}

      //! Constructor setting up a vector iniitialized with a constant value
      DynamicVector ( unsigned int size, const FieldType &s )
      : fields_( size, s )
      {}

      //! Copy constructor setting up a vector with the data of another one
      template< class T >
      DynamicVector ( const VectorInterface< T > &v )
      : fields_()
      {
        assign( v );
      }

      //! Copy constructor setting up a vector with the data of another one (of the same type)
      DynamicVector ( const ThisType &v )
      : fields_()
      {
        assign( v );
      }

      //! Assign another vector to this one
      template< class T >
      ThisType &operator= ( const VectorInterface< T > &v )
      {
        assign( v );
        return *this;
      }

      //! Assign another vector (of the same type) to this one
      ThisType &operator= ( const ThisType &v )
      {
        assign( v );
        return *this;
      }

      //! Initialize all fields of this vector with a scalar
      ThisType &operator= ( const FieldType &s )
      {
        assign( s );
        return *this;
      }

      const FieldType &operator[] ( unsigned int index ) const
      {
        return fields_[ index ];
      }

      FieldType &operator[] ( unsigned int index )
      {
        return fields_[ index ];
      }

      /** \copydoc Dune::Fem::VectorInterface::assign(const VectorInterface<T> &v) */
      template< class T >
      void assign ( const VectorInterface< T > &v )
      {
        fields_.assign( v.begin(), v.end() );
      }

      const FieldType *leakPointer () const
      {
        return fields_.data();
      }

      FieldType *leakPointer ()
      {
        return fields_.data();
      }

      void reserve ( unsigned int newSize )
      {
        fields_.reserve( newSize );
      }

      void resize ( unsigned int newSize )
      {
        fields_.resize( newSize );
      }

      void resize ( unsigned int newSize, const FieldType &defaultValue )
      {
        fields_.resize( newSize, defaultValue );
      }

      unsigned int size () const
      {
        return fields_.size();
      }
    };



    /** \class StaticVector
     *  \brief A vector using a std::array as storage
     *
     *  An implementation of VectorInterface using a std::array to provide the
     *  fields.
     */
    template< class FieldImp, int sz >
    class StaticVector
    : public VectorDefault< FieldImp, StaticVector< FieldImp, sz > >
    {
    public:
      //! field type of vector
      typedef FieldImp FieldType;

    private:
      typedef StaticVector< FieldImp, sz > ThisType;
      typedef VectorDefault< FieldImp, ThisType > BaseType;

    public:
      using BaseType :: assign;

    protected:
      std::array<FieldType,sz> fields_;

    public:
      //! Constructor setting up an uninitialized vector
      StaticVector () = default;

      //! Constructor setting up a vector initialized to a constant value
      explicit StaticVector ( const FieldType &s )
      {
        assign( s );
      }

      //! Copy constructor setting up a vector with the data of another one
      template< class T >
      StaticVector ( const VectorInterface< T > &v )
      {
        assign( v );
      }

      //! Copy constructor setting up a vector with the data of another one
      StaticVector ( const ThisType &v )
      {
        assign( v );
      }

      //! Assign another vector to this one
      template< class T >
      ThisType &operator= ( const VectorInterface< T > &v )
      {
        assign( v );
        return *this;
      }

      //! Assign another vector to this one
      ThisType &operator= ( const ThisType &v )
      {
        assign( v );
        return *this;
      }

      //! Initialize all fields of this vector with a scalar
      ThisType &operator= ( const FieldType &s )
      {
        assign( s );
        return *this;
      }

      const FieldType &operator[] ( unsigned int index ) const
      {
        assert( index < sz );
        return fields_[ index ];
      }

      FieldType &operator[] ( unsigned int index )
      {
        assert( index < sz );
        return fields_[ index ];
      }

      unsigned int size () const
      {
        return sz;
      }
    };



    // Capabilities
    // ------------
    namespace Capabilities
    {

      template< class Array >
      struct HasLeakPointer
      : public MetaBool< false >
      {};

      template< class Field, template< class > class Allocator >
      struct HasLeakPointer< std::vector< Field, Allocator<Field> > >
      : public MetaBool< true >
      {};

    }

  } // namespace Fem

} // namespace Dune

#include "vector_inline.hh"

//! @}

#endif // #ifndef DUNE_FEM_VECTOR_HH
