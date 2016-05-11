#ifndef DUNE_FEM_VECTOR_HH
#define DUNE_FEM_VECTOR_HH

#include <algorithm>
#include <cassert>
#include <type_traits>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/misc/bartonnackmaninterface.hh>


/*! @addtogroup VectorClasses
    @{
*/

namespace Dune
{
  namespace Fem
  {
    /** \class VectorInterface
     *  \ingroup Vector
     *  \brief Abstract vector interface
     */
    template< class VT >
    class VectorInterface
    : public BartonNackmanInterface< VectorInterface< VT >, typename VT::VectorType >
    {
      typedef VectorInterface< VT > ThisType;
      typedef BartonNackmanInterface< ThisType, typename VT::VectorType > BaseType;

    public:
      //! Type of the traits
      typedef VT Traits;

      //! Type of this interface
      typedef ThisType VectorInterfaceType;

      //! Type of the implementation (Barton-Nackman)
      typedef typename Traits::VectorType VectorType;

      //! Field type for the vector
      typedef typename Traits::FieldType FieldType;
      typedef FieldType value_type;

      //! Type of constant iterator
      typedef typename Traits::ConstIteratorType ConstIteratorType;
      typedef ConstIteratorType const_iterator;

      //! Type of iterator
      typedef typename Traits::IteratorType IteratorType;
      typedef IteratorType iterator;

      //! Type of unsigned integral type of indexing
      typedef unsigned int  size_type;

      //! Assign another vector to this one
      template< class T >
      VectorType& operator= ( const VectorInterface< T > &v )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( v ) );
        return asImp();
      }

      //! Assign another vector to this one
      VectorType& operator= ( const ThisType &v )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( v ) );
        return asImp();
      }

      //! Initialize all fields of this vector with a scalar
      VectorType &operator= ( const FieldType &s )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( s ) );
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

      //! Copy another vector to this one
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

      //! Initialize the vector to 0
      void clear ()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().clear() );
      }

      //! Obtain begin iterator
      ConstIteratorType begin () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      //! Obtain begin iterator
      IteratorType begin ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      //! Obtain end iterator
      ConstIteratorType end () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      //! Obtain end iterator
      IteratorType end ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      //! Obtain vector's size
      unsigned int size () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
        return asImp().size();
      }

    protected:
      using BaseType::asImp;
    };



    template< class Element, class Vector >
    class VectorDefaultIterator
      : public ForwardIteratorFacade< VectorDefaultIterator< Element, Vector >, Element >
    {
      typedef VectorDefaultIterator< Element, Vector > ThisType;

    public:
      typedef Element ElementType;
      typedef Vector VectorType;

      VectorDefaultIterator ( VectorType &vector, unsigned int index )
      : vector_( vector ), index_( index )
      {
        assert( index <= vector.size() );
      }

      VectorDefaultIterator( const ThisType &other ) = default;

      ThisType &operator= ( const ThisType &other )
      {
        assert( &(other.vector_) == &vector_ );
        index_ = other.index_;
      }

      ElementType& dereference() const
      {
        assert( index_ < vector_.size() );
        return vector_[ index_ ];
      }

      void increment()
      {
        assert( index_ < vector_.size() );
        ++index_;
      }

      bool equals( const ThisType &other ) const
      {
        assert( &(other.vector_) == &vector_ );
        return index_ == other.index_;
      }

      unsigned int index () const
      {
        return index_;
      }

    private:
      VectorType &vector_;
      unsigned int index_;
    };



    template< class Field, class Vector >
    struct VectorDefaultTraits
    {
      typedef Field FieldType;
      typedef Vector VectorType;

      typedef VectorDefaultIterator< FieldType, VectorType > IteratorType;
      typedef VectorDefaultIterator< const FieldType, const VectorType > ConstIteratorType;
    };



    template< class Vector >
    struct SupportsVectorInterface
    {
      typedef VectorInterface< typename Vector::Traits > VectorInterfaceType;
      static const bool v = std::is_convertible< Vector, VectorInterfaceType >::value;
    };



    /** \class VectorDefault
     *  \ingroup Vector
     *  \brief Default implementation of VectorInterface
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

      /** \copydoc Dune::Fem::VectorInterface::operator=(const FieldType &s) */
      VectorType &operator= ( const FieldType &s )
      {
        asImp().assign( s );
        return asImp();
      }

      /** \copydoc Dune::Fem::VectorInterface::operator=(const VectorInterface<T> &v) */
      template< class T >
      ThisType &operator= ( const VectorInterface< T > &v )
      {
        asImp().assign( v );
        return asImp();
      }

      /** \copydoc Dune::Fem::VectorInterface::operator+=(const VectorInterface<T> &v) */
      template< class T >
      VectorType &operator+= ( const VectorInterface< T > &v )
      {
        const unsigned int size = this->size();
        assert( size == v.size() );
        for( unsigned int i = 0; i < size; ++i )
          (*this)[ i ] += v[ i ];
        return asImp();
      }

      /** \copydoc Dune::Fem::VectorInterface::operator-=(const VectorInterface<T> &v) */
      template< class T >
      VectorType &operator-= ( const VectorInterface< T > &v )
      {
        const unsigned int size = this->size();
        assert( size == v.size() );
        for( unsigned int i = 0; i < size; ++i )
          (*this)[ i ] -= v[ i ];
        return asImp();
      }

      /** \copydoc Dune::Fem::VectorInterface::operator*=(const FieldType &s) */
      VectorType &operator*= ( const FieldType &s )
      {
        for( auto& entry : *this )
          entry *= s;
        return asImp();
      }

      /** \copydoc Dune::Fem::VectorInterface::addScaled(const FieldType &s,const VectorInterface<T> &v) */
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

      /** \copydoc Dune::Fem::VectorInterface::assign(const FieldType &s) */
      void assign ( const FieldType &s )
      {
        std::fill( asImp().begin(), asImp().end(), s );
      }

      /** \copydoc Dune::Fem::VectorInterface::clear() */
      void clear ()
      {
        asImp().assign( 0 );
      }

      /** \copydoc Dune::Fem::VectorInterface::begin() */
      ConstIteratorType begin () const
      {
        return ConstIteratorType( asImp(), 0 );
      }

      /** \copydoc Dune::Fem::VectorInterface::begin() */
      IteratorType begin ()
      {
        return IteratorType( asImp(), 0 );
      }

      /** \copydoc Dune::Fem::VectorInterface::end() */
      ConstIteratorType end () const
      {
        return ConstIteratorType( asImp(), size() );
      }

      /** \copydoc Dune::Fem::VectorInterface::end() */
      IteratorType end ()
      {
        return IteratorType( asImp(), size() );
      }

      using BaseType :: size;

    protected:
      using BaseType :: asImp;
    };

  } // namespace Fem
} // namespace Dune

#include "vector_inline.hh"

//! @}

#endif // #ifndef DUNE_FEM_VECTOR_HH
