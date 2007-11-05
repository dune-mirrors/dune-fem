#ifndef DUNE_FEM_VECTOR_HH
#define DUNE_FEM_VECTOR_HH

#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/storage/arrayallocator.hh>
#include <dune/fem/storage/array.hh>

namespace Dune
{

  template< class VectorTraits >
  struct VectorInterfaceArrayTraits
  {
    typedef typename VectorTraits :: VectorType ArrayType;

    typedef typename VectorTraits :: FieldType ElementType;

    typedef typename VectorTraits :: ConstIteratorType ConstIteratorType;
    
    typedef typename VectorTraits :: IteratorType IteratorType;
  };



  //! An abstract vector interface
  template< class TraitsImp >
  class VectorInterface
  : public ArrayInterface< VectorInterfaceArrayTraits< TraitsImp > >
  {
  public:
    typedef TraitsImp Traits;

  private:
    typedef VectorInterface< Traits > ThisType;
    typedef ArrayInterface< VectorInterfaceArrayTraits< Traits > > BaseType;

    template< class >
    friend class VectorInterface;

  public:
    //! type of this interface
    typedef ThisType VectorInterfaceType;

    //! type of the implementation (Barton-Nackman)
    typedef typename Traits :: VectorType VectorType;

    //! field type for the vector
    typedef typename Traits :: FieldType FieldType;

    //! type of constant iterator
    typedef typename Traits :: ConstIteratorType ConstIteratorType;
    
    //! type of iterator
    typedef typename Traits :: IteratorType IteratorType;

  protected:
    using BaseType :: asImp;

  public:
    //! Assign another vector to this one
    template< class T >
    inline VectorType& operator= ( const VectorInterface< T > &v )
    {
      return asImp().assign( v );
    }
    
    //! Assign another vector to this one
    inline VectorType& operator= ( const ThisType &v )
    {
      return asImp().assign( v );
    }

    //! Initialize all fields of this vector with a scalar
    inline VectorType &operator= ( const FieldType s )
    {
      return asImp().assign( s );
    }

    //! Returns a const reference to the field indexed by index
    inline const FieldType &operator[] ( unsigned int index ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
    }

    //! Returns a reference to the field indexed by index
    inline FieldType &operator[] ( unsigned int index )
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
    }

    //! Add another vector to this one
    template< class T >
    inline VectorType &operator+= ( const VectorInterface< T > &v )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator+=( v.asImp() ) );
      return asImp();
    }

    //! Subtract another vector from this one
    template< class T >
    inline VectorType &operator-= ( const VectorInterface< T > &v )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator-=( v.asImp() ) );
      return asImp();
    }

    //! Multiply this vector by a scalar
    inline VectorType &operator*= ( const FieldType s )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator*=( s ) );
      return asImp();
    }
            
    //! Add a multiple of another vector to this one
    template< class T >
    inline VectorType &addScaled ( const FieldType s, const VectorInterface< T > &v )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().add( s, v.asImp() ) );
      return asImp();
    }
    
    //! Assign another vector to this one
    template< class T >
    inline VectorType &assign ( const VectorInterface< T > &v )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( v.asImp() ) );
      return asImp();
    }

    //! Initialize all fields of this vector with a scalar
    inline VectorType &assign ( const FieldType s )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( s ) );
      return asImp();
    }

    //! obtain begin iterator
    inline ConstIteratorType begin () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
      return asImp().begin();
    }

    //! obtain begin iterator
    inline IteratorType begin ()
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
      return asImp().begin();
    }

    //! obtain end iterator
    inline ConstIteratorType end () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
      return asImp().end();
    }

    //! obtain end iterator
    inline IteratorType end ()
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
      return asImp().end();
    }
    
    //! Returns the vector's size
    inline unsigned int size () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
      return asImp().size();
    }
  };



  template< class VectorType >
  struct CheckVectorInterface
  {
    typedef VectorInterface< typename VectorType :: Traits > VectorInterfaceType;
    
    typedef CompileTimeChecker< Conversion< VectorType, VectorInterfaceType > :: exists >
      CheckerType;
  };



  template< class Type1, class Type2 >
  struct ExtractCommonFieldType
  {
    typedef typename Type1 :: FieldType FieldType;

    typedef CompileTimeChecker
      < Conversion< FieldType, typename Type2 :: FieldType > :: sameType >
      __FieldType_Must_Be_Identical__;
  };



  template< class FieldImp, class VectorImp >
  struct VectorDefaultTraits
  {
    typedef FieldImp FieldType;

    typedef VectorImp VectorType;

    typedef ArrayDefaultIterator< FieldType, VectorType > IteratorType;

    typedef ArrayDefaultIterator< const FieldType, const VectorType > ConstIteratorType;
  };



  //! Default implementation of VectorInterface
  template< class FieldImp, class VectorImp >
  class VectorDefault
  : public VectorInterface< VectorDefaultTraits< FieldImp, VectorImp > >
  {
  public:
    typedef FieldImp FieldType;

    typedef VectorImp VectorType;

    typedef VectorDefaultTraits< FieldType, VectorType > Traits;

  private:
    typedef VectorDefault< FieldType, VectorType > ThisType;
    typedef VectorInterface< Traits > BaseType;

  protected:
    using BaseType :: asImp;

  public:
    using BaseType :: size;

  public:
    typedef typename BaseType :: VectorInterfaceType VectorInterfaceType;

    typedef typename Traits :: ConstIteratorType ConstIteratorType;
    typedef typename Traits :: IteratorType IteratorType;

  public:
    //! Add another vector to this one
    template< class T >
    inline VectorType &operator+= ( const VectorInterface< T > &v )
    {
      const unsigned int size = this->size();
      assert( size == v.size() );
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] += v[ i ];
      return asImp();
    }
   
    //! Subtract another vector from this one
    template< class T >
    inline VectorType &operator-= ( const VectorInterface< T > &v )
    {
      const unsigned int size = this->size();
      assert( size == v.size() );
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] -= v[ i ];
      return asImp();
    }

    //! Multiply this vector by a scalar
    inline VectorType &operator*= ( const FieldType s )
    {
      const unsigned int size = this->size();
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] *= s;
      return asImp();
    }
    
    //! Add a multiple of another vector to this one
    template< class T >
    inline VectorType &addScaled ( const FieldType s,
                                   const VectorInterface< T > &v )
    {
      const unsigned int size = this->size();
      assert( size == v.size() );
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] += s * v[ i ];
      return asImp();
    }

    //! Assign another vector to this one
    template< class T >
    inline VectorType &assign ( const VectorInterface< T > &v )
    {
      const unsigned int size = this->size();
      assert( size == v.size() );
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] = v[ i ];
      return asImp();
    }
    
    //! Initialize all fields of this vector with a scalar
    inline VectorType &assign ( const FieldType s )
    {
      const unsigned int size = this->size();
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] = s;
      return asImp();
    }

    //! obtain begin iterator
    inline ConstIteratorType begin () const
    {
      return ConstIteratorType( asImp(), 0 );
    }

    //! obtain begin iterator
    inline IteratorType begin ()
    {
      return IteratorType( asImp(), 0 );
    }

    //! obtain end iterator
    inline ConstIteratorType end () const
    {
      return ConstIteratorType( asImp(), size() );
    }

    //! obtain end iterator
    inline IteratorType end ()
    {
      return IteratorType( asImp(), size() );
    }
  };



  template< class FieldVectorImp >
  class FieldVectorWrapper;



  template< class FieldImp, int sz >
  class FieldVectorWrapper< FieldVector< FieldImp, sz > >
  : public VectorDefault< FieldImp, FieldVectorWrapper< FieldVector < FieldImp, sz > > >
  {
  public:
    typedef FieldImp FieldType;

    typedef FieldVector< FieldType, sz > FieldVectorType;

  private:
    typedef FieldVectorWrapper< FieldVectorType > ThisType;
    typedef VectorDefault< FieldType, ThisType > BaseType;

  public:
    using BaseType :: operator+=;
    using BaseType :: operator-=;
    using BaseType :: addScaled;
    using BaseType :: assign;

  protected:
    FieldVectorType fieldVector_;

  public:
    inline FieldVectorWrapper ()
    : fieldVector_()
    {
    }

    inline explicit FieldVectorWrapper ( const FieldType s )
    : fieldVector_( s )
    {
    }

    inline explicit FieldVectorWrapper ( const FieldVectorType &v )
    : fieldVector_( v )
    {
    }

    template< class T >
    inline FieldVectorWrapper ( const VectorInterface< T > &v )
    : fieldVector_()
    {
      assign( v );
    }
    
    inline FieldVectorWrapper ( ThisType &other )
    : fieldVector_( other.fieldVector_ )
    {
    }

  public:
    inline operator const FieldVectorType& () const
    {
      return fieldVector_;
    }

    inline operator FieldVectorType& ()
    {
      return fieldVector_;
    }

    template< class T >
    inline ThisType &operator= ( const VectorInterface< T > &v )
    {
      return assign( v );
    }

    inline ThisType &operator= ( const ThisType &v )
    {
      return assign( v );
    }

    inline ThisType &operator= ( const FieldType s )
    {
      return assign( s );
    }
    
    inline const FieldType &operator[] ( unsigned int index ) const
    {
      return fieldVector_[ index ];
    }

    inline FieldType &operator[] ( unsigned int index )
    {
      return fieldVector_[ index ];
    }

    inline ThisType &operator+= ( const ThisType &v )
    {
      fieldVector_ += v.fieldVector_;
      return *this;
    }

    inline ThisType &operator+= ( const FieldVectorType &v )
    {
      fieldVector_ += v;
      return *this;
    }

    inline ThisType &operator-= ( const ThisType &v )
    {
      fieldVector_ += v.fieldVector_;
      return *this;
    }

    inline ThisType &operator-= ( const FieldVectorType &v )
    {
      fieldVector_ -= v;
      return *this;
    }

    inline ThisType &operator*= ( const FieldType s )
    {
      fieldVector_ *= s;
      return *this;
    }

    inline ThisType &addScaled ( const FieldType s, const ThisType &other )
    {
      fieldVector_.axpy( s, other.fieldVector_ );
      return *this;
    }
    
    inline ThisType &assign ( const ThisType &other )
    {
      fieldVector_ = other.fieldVector_;
      return *this;
    }

    inline ThisType &assign ( const FieldType s )
    {
      fieldVector_ = s;
      return *this;
    }

    inline unsigned int size () const
    {
      return FieldVectorType :: size;
    }
  };



  //! An implementation of VectorInterface wrapping a standard C++ array
  template< class FieldImp >
  class ArrayWrapperVector
  : public VectorDefault< FieldImp, ArrayWrapperVector< FieldImp > >
  {
  public:
    //! field type of vector
    typedef FieldImp FieldType;

  private:
    typedef ArrayWrapperVector< FieldType > ThisType;
    typedef VectorDefault< FieldType, ThisType > BaseType;

  public:
    using BaseType :: assign;

  protected:
    const unsigned int size_;
    FieldType *const fields_;

  public:
    //! Constructor setting up the vector (without initializing the fields)
    inline ArrayWrapperVector ( const unsigned int size,
                                FieldType *const fields )
    : size_( size ),
      fields_( fields )
    {
    }

    //! Constructor setting up the vector and initializing the fields to a constant value
    inline ArrayWrapperVector ( const unsigned int size,
                                FieldType *const fields,
                                const FieldType s )
    : size_( size ),
      fields_( fields )
    {
      assign( s );
    }

    //! Copy constructor setting up a vector with the data of another one
    template< class T >
    inline ArrayWrapperVector ( const unsigned int size,
                                FieldType *const fields,
                                const VectorInterface< T > &v )
    : size_( size ),
      fields_( fields )
    {
      assign( v );
    }

    //! Assign another vector to this one
    template< class T >
    inline ThisType &operator= ( const VectorInterface< T > &v )
    {
      return assign( v );
    }

    //! Assign another vector to this one
    inline ThisType &operator= ( const ThisType &v )
    {
      return assign( v );
    }

    //! Initialize all fields of this vector with a scalar
    inline ThisType &operator= ( const FieldType s )
    {
      return assign( s );
    }

    inline const FieldType &operator[] ( unsigned int index ) const
    {
      assert( index < size_ );
      return fields_[ index ];
    }

    inline FieldType &operator[] ( unsigned int index )
    {
      assert( index < size_ );
      return fields_[ index ];
    }

    inline unsigned int size () const
    {
      return size_;
    }
  };



  /** \class DynamicVector
   *  \brief A vector using a DynamicArray as storage
   * 
   *  An implementation of VectorInterface using a DynamicArray to provide the
   *  fields.
   */
  template< class FieldImp,
            template< class > class ArrayAllocatorImp = DefaultArrayAllocator >
  class DynamicVector
  : public VectorDefault< FieldImp, DynamicVector< FieldImp, ArrayAllocatorImp > >
  {
  public:
    //! field type of the vector
    typedef FieldImp FieldType;

  private:
    typedef DynamicVector< FieldType, ArrayAllocatorImp > ThisType;
    typedef VectorDefault< FieldType, ThisType > BaseType;

  public:
    using BaseType :: assign;
    
  protected:
    DynamicArray< FieldType, ArrayAllocatorImp > fields_;

  public:
    //! Constructor setting up a vector of a specified size
    inline explicit DynamicVector ( unsigned int size = 0 )
    : fields_( size )
    {
    }

    //! Constructor setting up a vector iniitialized with a constant value
    inline DynamicVector ( unsigned int size,
                           const FieldType s )
    : fields_( size )
    {
      assign( s );
    }

    //! Copy constructor setting up a vector with the data of another one
    template< class T >
    inline DynamicVector ( const VectorInterface< T > &v )
    : fields_()
    {
      assign( v );
    }

    //! Copy constructor setting up a vector with the data of another one (of the same type)
    inline DynamicVector ( const ThisType &v )
    : fields_()
    {
      assign( v );
    }

    //! Assign another vector to this one
    template< class T >
    inline ThisType &operator= ( const VectorInterface< T > &v )
    {
      return assign( v );
    }

    //! Assign another vector (of the same type) to this one
    inline ThisType &operator= ( const ThisType &v )
    {
      return assign( v );
    }

    //! Initialize all fields of this vector with a scalar
    inline ThisType &operator= ( const FieldType s )
    {
      return assign( s );
    }

    inline const FieldType &operator[] ( unsigned int index ) const
    {
      return fields_[ index ];
    }
    
    inline FieldType &operator[] ( unsigned int index )
    {
      return fields_[ index ];
    }
    
    /** \copydoc Dune::VectorInterface::assign(const VectorInterrace<T> &v) */
    template< class T >
    inline ThisType &assign ( const VectorInterface< T > &v )
    {
      fields_.assign( v );
      return *this;
    }

    inline void resize ( unsigned int newSize )
    {
      fields_.resize( newSize );
    }

    inline void resize ( unsigned int newSize,
                         const FieldType defaultValue )
    {
      fields_.resize( newSize, defaultValue );
    }

    inline unsigned int size () const
    {
      return fields_.size();
    }
  };



  /*! \class StaticVector
   *  \brief implementation of VectorInterface using a C++ array embedded info
   *         the class to provide the fields
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
    FieldType fields_[ sz ]; //!< The actual vector fields

  public:
    //! Constructor setting up an uninitialized vector
    inline StaticVector ()
    {
    }

    //! Constructor setting up a vector initialized to a constant value
    inline explicit StaticVector ( const FieldType s )
    {
      assign( s );
    }

    //! Copy constructor setting up a vector with the data of another one
    template< class T >
    inline StaticVector ( const VectorInterface< T > &v )
    {
      assign( v );
    }
    
    //! Copy constructor setting up a vector with the data of another one
    inline StaticVector ( const ThisType &v )
    {
      assign( v );
    }

    //! Assign another vector to this one
    template< class T >
    inline ThisType &operator= ( const VectorInterface< T > &v )
    {
      return assign( v );
    }

    //! Assign another vector to this one
    inline ThisType &operator= ( const ThisType &v )
    {
      return assign( v );
    }

    //! Initialize all fields of this vector with a scalar
    inline ThisType &operator= ( const FieldType s )
    {
      return assign( s );
    }

    inline const FieldType &operator[] ( unsigned int index ) const
    {
      assert( index < sz );
      return fields_[ index ];
    }

    inline FieldType &operator[] ( unsigned int index )
    {
      assert( index < sz );
      return fields_[ index ];
    }

    inline unsigned int size () const
    {
      return sz;
    }
  };


  
  template< class Vector1Type, class Vector2Type >
  class CombinedVector
  : public VectorDefault< typename Vector1Type :: FieldType,
                          CombinedVector< Vector1Type, Vector2Type > >
  {
  public:
    typedef typename Vector1Type :: FieldType FieldType;

  private:
    typedef CombinedVector< Vector1Type, Vector2Type > ThisType;
    typedef VectorDefault< typename Vector1Type :: FieldType, ThisType > BaseType;

    typedef CheckVectorInterface< Vector1Type > __CheckVector1Type__;
    typedef CheckVectorInterface< Vector2Type > __CheckVector2Type__;

    typedef CompileTimeChecker
      < Conversion< FieldType, typename Vector2Type :: FieldType > :: sameType >
      __CheckFieldType__;

  protected:
      Vector1Type &vector1_;
      Vector2Type &vector2_;

  public:
    inline CombinedVector( Vector1Type &v1, Vector2Type &v2 )
    : vector1_( v1 ),
      vector2_( v2 )
    {
    }
      
    inline const FieldType &operator[] ( unsigned int index ) const
    {
      const int index2 = index - vector1_.size();
      if( index2 < 0 )
        return vector1_[ index ];
      else
        return vector2_[ index2 ];
    }

    inline FieldType &operator[] ( unsigned int index )
    {
      const int index2 = index - vector1_.size();
      if( index2 < 0 )
        return vector1_[ index ];
      else
        return vector2_[ index2 ];
    }

    inline const unsigned int size() const
    {
      return vector1_.size() + vector2_.size();
    }
  };

}

#include "vector_inline.hh"

#endif
