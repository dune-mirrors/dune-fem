#ifndef DUNE_FEM_VECTOR_HH
#define DUNE_FEM_VECTOR_HH

#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/misc/combineinterface.hh>

#include <dune/fem/storage/arrayallocator.hh>
#include <dune/fem/storage/array.hh>
#include <dune/fem/io/streams/streams.hh>

/*! @addtogroup VectorClasses
    @{
*/

namespace Dune
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

  public:
    //! Assign another vector to this one
    template< class T >
    VectorType& operator= ( const VectorInterface< T > &v );
    
    //! Assign another vector to this one
    VectorType& operator= ( const ThisType &v );

    //! Initialize all fields of this vector with a scalar
    VectorType &operator= ( const FieldType s );

    //! Returns a const reference to the field indexed by index
    const FieldType &operator[] ( unsigned int index ) const;

    //! Returns a reference to the field indexed by index
    FieldType &operator[] ( unsigned int index );

    //! Add another vector to this one
    template< class T >
    VectorType &operator+= ( const VectorInterface< T > &v );

    //! Subtract another vector from this one
    template< class T >
    VectorType &operator-= ( const VectorInterface< T > &v );

    //! Multiply this vector by a scalar
    VectorType &operator*= ( const FieldType s );
            
    //! Add a multiple of another vector to this one
    template< class T >
    VectorType &addScaled ( const FieldType s,
                            const VectorInterface< T > &v );
    
    /** \brief copy another vector to this one
     *
     *  Copies the data from another vector to this one. Both vectors must be of
     *  the same size.
     *
     *  \param[in]  v  vector to copy
     */
    template< class T >
    void assign ( const VectorInterface< T > &v );

    //! Initialize all fields of this vector with a scalar
    void assign ( const FieldType s );

    /** \brief initialize the vector to 0 */
    void clear ();

    //! obtain begin iterator
    ConstIteratorType begin () const;

    //! obtain begin iterator
    IteratorType begin ();

    //! obtain end iterator
    ConstIteratorType end () const;

    //! obtain end iterator
    IteratorType end ();
    
    //! Returns the vector's size
    unsigned int size () const;

  protected:
    using BaseType::asImp;
  };


  template< class Vector >
  struct SupportsVectorInterface
  {
    typedef VectorInterface< typename Vector::Traits > VectorInterfaceType;
    static const bool v = Conversion< Vector, VectorInterfaceType >::exists;
  };


  template< class V, class W >
  struct ExtractCommonFieldType
  {
    typedef typename V::FieldType FieldType;

  private:
    dune_static_assert( (Conversion< FieldType, typename W::FieldType >::sameType),
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

  public:
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
    VectorType &operator*= ( const FieldType s )
    {
      const unsigned int size = this->size();
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] *= s;
      return asImp();
    }
    
    //! Add a multiple of another vector to this one
    template< class T >
    VectorType &addScaled ( const FieldType s, const VectorInterface< T > &v )
    {
      const unsigned int size = this->size();
      assert( size == v.size() );
      for( unsigned int i = 0; i < size; ++i )
        (*this)[ i ] += s * v[ i ];
      return asImp();
    }

    /** \copydoc Dune::VectorInterface::assign(const VectorInterface<T> &v) */
    template< class T >
    void assign ( const VectorInterface< T > &v )
    {
      const unsigned int size = this->size();
      assert( size == v.size() );
      for( unsigned int i = 0; i < size; ++i )
        asImp()[ i ] = v[ i ];
    }
    
    //! Initialize all fields of this vector with a scalar
    void assign ( const FieldType s )
    {
      const unsigned int size = this->size();
      for( unsigned int i = 0; i < size; ++i )
        asImp()[ i ] = s;
    }

    /** \copydoc Dune::VectorInterface::clear() */
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



  // FieldVectorAdapter
  // ------------------

  template< class FieldVectorImp >
  class FieldVectorAdapter;

  template< class Field, int sz >
  class FieldVectorAdapter< FieldVector< Field, sz > >
  : public VectorDefault< Field, FieldVectorAdapter< FieldVector< Field, sz > > >
  {
    typedef FieldVectorAdapter< FieldVector< Field, sz > > ThisType;
    typedef VectorDefault< Field, ThisType > BaseType;

  public:
    typedef Field FieldType;

    typedef FieldVector< FieldType, sz > FieldVectorType;

    using BaseType :: operator+=;
    using BaseType :: operator-=;
    using BaseType :: addScaled;
    using BaseType :: assign;

  protected:
    FieldVectorType fieldVector_;

  public:
    FieldVectorAdapter ()
    : fieldVector_()
    {}

    explicit FieldVectorAdapter ( const FieldType &s )
    : fieldVector_( s )
    {}

    explicit FieldVectorAdapter ( const FieldVectorType &v )
    : fieldVector_( v )
    {}

    template< class T >
    FieldVectorAdapter ( const VectorInterface< T > &v )
    : fieldVector_()
    {
      assign( v );
    }
    
    FieldVectorAdapter ( const ThisType &other )
    : fieldVector_( other.fieldVector_ )
    {}
 
  public:
    operator const FieldVectorType & () const
    {
      return fieldVector_;
    }

    operator FieldVectorType & ()
    {
      return fieldVector_;
    }

    template< class T >
    ThisType &operator= ( const VectorInterface< T > &v )
    {
      assign( v );
      return *this;
    }

    ThisType &operator= ( const ThisType &v )
    {
      assign( v );
      return *this;
    }

    ThisType &operator= ( const FieldType &s )
    {
      return assign( s );
    }
    
    const FieldType &operator[] ( unsigned int index ) const
    {
      return fieldVector_[ index ];
    }

    FieldType &operator[] ( unsigned int index )
    {
      return fieldVector_[ index ];
    }

    ThisType &operator+= ( const ThisType &v )
    {
      fieldVector_ += v.fieldVector_;
      return *this;
    }

    ThisType &operator+= ( const FieldVectorType &v )
    {
      fieldVector_ += v;
      return *this;
    }

    ThisType &operator-= ( const ThisType &v )
    {
      fieldVector_ += v.fieldVector_;
      return *this;
    }

    ThisType &operator-= ( const FieldVectorType &v )
    {
      fieldVector_ -= v;
      return *this;
    }

    ThisType &operator*= ( const FieldType &s )
    {
      fieldVector_ *= s;
      return *this;
    }

    ThisType &addScaled ( const FieldType &s, const ThisType &other )
    {
      fieldVector_.axpy( s, other.fieldVector_ );
      return *this;
    }
    
    void assign ( const ThisType &other )
    {
      fieldVector_ = other.fieldVector_;
    }

    void assign ( const FieldType &s )
    {
      fieldVector_ = s;
    }

    unsigned int size () const
    {
      return FieldVectorType :: size;
    }

    static const ThisType &adapt ( const FieldVectorType &v )
    {
      return reinterpret_cast< const ThisType & >( v );
    }

    static ThisType &adapt ( FieldVectorType &v )
    {
      return reinterpret_cast< ThisType & >( v );
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
      assign( v );
      return *this;
    }

    //! Assign another vector to this one
    inline ThisType &operator= ( const ThisType &v )
    {
      assign( v );
      return *this;
    }

    //! Initialize all fields of this vector with a scalar
    inline ThisType &operator= ( const FieldType s )
    {
      assign( s );
      return *this;
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
  template< class Field,
            template< class > class ArrayAllocator = DefaultArrayAllocator >
  class DynamicVector
  : public VectorDefault< Field, DynamicVector< Field, ArrayAllocator > >
  {
    typedef DynamicVector< Field, ArrayAllocator > ThisType;
    typedef VectorDefault< Field, ThisType > BaseType;

  public:
    //! field type of the vector
    typedef Field FieldType;

    using BaseType :: assign;
    
  protected:
    DynamicArray< FieldType, ArrayAllocator > fields_;

  public:
    //! Constructor setting up a vector of a specified size
    inline explicit DynamicVector ( unsigned int size = 0 )
    : fields_( size )
    {}

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
      assign( v );
      return *this;
    }

    //! Assign another vector (of the same type) to this one
    inline ThisType &operator= ( const ThisType &v )
    {
      assign( v );
      return *this;
    }

    //! Initialize all fields of this vector with a scalar
    inline ThisType &operator= ( const FieldType s )
    {
      assign( s );
      return *this;
    }

    inline const FieldType &operator[] ( unsigned int index ) const
    {
      return fields_[ index ];
    }
    
    inline FieldType &operator[] ( unsigned int index )
    {
      return fields_[ index ];
    }
    
    /** \copydoc Dune::VectorInterface::assign(const VectorInterface<T> &v) */
    template< class T >
    inline void assign ( const VectorInterface< T > &v )
    {
      fields_.assign( v );
    }

    inline const FieldType *leakPointer () const
    {
      return fields_.leakPointer();
    }

    inline FieldType *leakPointer ()
    {
      return fields_.leakPointer();
    }
    
    inline void reserve ( unsigned int newSize )
    {
      fields_.reserve( newSize );
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
      assign( v );
      return *this;
    }

    //! Assign another vector to this one
    inline ThisType &operator= ( const ThisType &v )
    {
      assign( v );
      return *this;
    }

    //! Initialize all fields of this vector with a scalar
    inline ThisType &operator= ( const FieldType s )
    {
      assign( s );
      return *this;
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
  : public VectorDefault< typename ExtractCommonFieldType< Vector1Type, Vector2Type >::FieldType,
                          CombinedVector< Vector1Type, Vector2Type > >
  {
    typedef CombinedVector< Vector1Type, Vector2Type > ThisType;
    typedef VectorDefault< typename ExtractCommonFieldType< Vector1Type, Vector2Type >::FieldType, ThisType > BaseType;

    dune_static_assert( SupportsVectorInterface< Vector1Type >::v, "CombinedVector only works on vectors." );
    dune_static_assert( SupportsVectorInterface< Vector2Type >::v, "CombinedVector only works on vectors." );

  public:
    typedef typename ExtractCommonFieldType< Vector1Type, Vector2Type >::FieldType FieldType;
    
  public:
    CombinedVector( Vector1Type &v1, Vector2Type &v2 )
    : vector1_( v1 ),
      vector2_( v2 )
    {}
      
    const FieldType &operator[] ( unsigned int index ) const
    {
      const int index2 = index - vector1_.size();
      if( index2 < 0 )
        return vector1_[ index ];
      else
        return vector2_[ index2 ];
    }

    FieldType &operator[] ( unsigned int index )
    {
      const int index2 = index - vector1_.size();
      if( index2 < 0 )
        return vector1_[ index ];
      else
        return vector2_[ index2 ];
    }

    const unsigned int size() const
    {
      return vector1_.size() + vector2_.size();
    }

  protected:
    Vector1Type &vector1_;
    Vector2Type &vector2_;
  };



  /** \class PairOfVector
   *  \brief Combines two classes implementing the vector interface to
   *         produce a long vector.

   *  This class can be used together with Dune::CombinedInterface
   *  to generate a large vector by combining smaller ones. 
   *  The class defines the Dune::VectorInterface
   *
   *  Example of usage:
   *  \code
      Dune::StaticVector<double,10> v1;
      Dune::DynamicVector<double> v2(2);
      Dune::FieldVectorWrapper<Dune::FieldVector<double,5> > v3;
      typedef Dune::CombineInterface<PairOfVectors,
            Dune::StaticVector<double,10>&,
            Dune::DynamicVector<double>&, 
            Dune::FieldVectorWrapper<Dune::FieldVector<double,5> >&
            VectorOfVectorType;
      VectorOfVectorType v(v1,v2,v3);
   *  \endcode
   */
  template< class RVector1Type, class RVector2Type >
  class PairOfVectors
  : public VectorDefault< typename TypeTraits< RVector1Type >::ReferredType::FieldType,
                          PairOfVectors<RVector1Type,RVector2Type> >,
    public PairOfInterfaces<RVector1Type,RVector2Type>
  {
    typedef PairOfVectors< RVector1Type, RVector2Type > ThisType;
    typedef VectorDefault< typename TypeTraits< RVector1Type >::ReferredType::FieldType, ThisType > BaseType;

    typedef typename PairOfInterfaces< RVector1Type, RVector2Type >::T1Type Vector1Type;
    typedef typename PairOfInterfaces< RVector1Type, RVector2Type >::T2Type Vector2Type;

    dune_static_assert( SupportsVectorInterface< Vector1Type >::v, "PairOfVectors only works on vectors." );
    dune_static_assert( SupportsVectorInterface< Vector2Type >::v, "PairOfVectors only works on vectors." );

  public:
    typedef typename ExtractCommonFieldType< Vector1Type, Vector2Type >::FieldType FieldType;

  private:
    PairOfVectors ();
    PairOfVectors &operator= ( const ThisType & );

  public:
    PairOfVectors ( RVector1Type v1, RVector2Type v2 )
    : Dune::PairOfInterfaces< RVector1Type, RVector2Type >( v1, v2 )
    {}
      
    const FieldType &operator[] ( unsigned int index ) const
    {
      const int index2 = index - this->first().size();
      if( index2 < 0 )
        return this->first()[ index ];
      else
        return this->second()[ index2 ];
    }

    FieldType &operator[] ( unsigned int index )
    {
      const int index2 = index - this->first().size();
      if( index2 < 0 )
        return this->first()[ index ];
      else
        return this->second()[ index2 ];
    }

    const unsigned int size () const
    {
      return this->first().size() + this->second().size();
    }
  };

  

  // Capabilities
  // ------------

  namespace Capabilities
  {

    template< class Field, template< class > class ArrayAllocator >
    struct HasLeakPointer< DynamicVector< Field, ArrayAllocator > >
    : public MetaBool< true >
    {};
    
  }
  
}

#include "vector_inline.hh"

//! @}

#endif // #ifndef DUNE_FEM_VECTOR_HH
