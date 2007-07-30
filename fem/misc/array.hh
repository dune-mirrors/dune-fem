#ifndef DUNE_FEM_ARRAY_HH
#define DUNE_FEM_ARRAY_HH

#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  /*! \class ArrayInterface
   *  \brief abstract array interface
   */
  template< class TraitsImp >
  class ArrayInterface
  {
  public:
    typedef TraitsImp TraitsType;

  private:
    typedef ArrayInterface< TraitsType > ThisType;

  public:
    //! type of this interface
    typedef ThisType ArrayInterfaceType;
    
    //! type of the implementation (Barton-Nackman) 
    typedef typename TraitsType :: ArrayType ArrayType;

  public:
    //! type of the array elements
    typedef typename TraitsType :: ElementType ElementType;

    //! type of constant iterator
    typedef typename TraitsType :: ConstIteratorType ConstIteratorType;

    //! type of iterator
    typedef typename TraitsType :: IteratorType IteratorType;

  public:
    inline ArrayInterface ()
    {
      typedef CompileTimeChecker< Conversion< ArrayType, ThisType > :: exists >
        __Array_Implementation_Must_Be_Derived_From_Interface__;
    }
    
    //! fill the array with copies of an element
    inline ArrayType &operator= ( const ElementType &element )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator=( element ) );
      return asImp();
    }
    
    //! access an array element
    inline const ElementType& operator[] ( unsigned int index ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
    }

    //! access an array element
    inline ElementType& operator[] ( unsigned int index )
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
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

    //! return the size of the array
    inline unsigned int size () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
      return asImp().size();
    }
    
  protected:
    inline const ArrayType& asImp () const
    {
      return static_cast< const ArrayType& >( *this );
    }

    inline ArrayType& asImp ()
    {
      return static_cast< ArrayType& >( *this );
    }
  };



  template< class ElementImp, class ArrayImp >
  class ArrayDefaultIterator
  {
  public:
    typedef ElementImp ElementType;

    typedef ArrayImp ArrayType;

  private:
    typedef ArrayDefaultIterator< ElementType, ArrayType > ThisType;

  protected:
    ArrayType &array_;
    unsigned int index_;

  public:
    inline ArrayDefaultIterator ( ArrayType &array,
                                  unsigned int index )
    : array_( array ),
      index_( index )
    {
      assert( index <= array.size() );
    }

    inline ArrayDefaultIterator( const ThisType &other )
    : array_( other.array_ ),
      index_( other.index_ )
    {
    }

    inline ThisType &operator= ( const ThisType &other )
    {
      assert( &(other.array_) == &array_ );
      index_ = other.index_;
    }

    inline ElementType &operator* ()
    {
      assert( index_ < array_.size() );
      return array_[ index_ ];
    }

    inline ThisType &operator++ ()
    {
      assert( index_ < array_.size() );
      ++index_;
      return *this;
    }

    inline bool operator== ( const ThisType &other )
    {
      assert( &(other.array_) == &array_ );
      return index_ == other.index_;
    }

    inline bool operator!= ( const ThisType &other )
    {
      return !(*this == other);
    }
  };



  template< class ElementImp, class ArrayImp >
  class ArrayDefaultTraits
  {
  public:
    typedef ElementImp ElementType;

    typedef ArrayImp ArrayType;

  private:
    typedef ArrayDefaultTraits< ElementType, ArrayType > ThisType;

  public:
    typedef ArrayDefaultIterator< ElementType, ArrayType > IteratorType;
    
    typedef ArrayDefaultIterator< const ElementType, const ArrayType > ConstIteratorType;
  };


  
  template< class ElementImp, class ArrayImp >
  class ArrayDefault
  : public ArrayInterface< ArrayDefaultTraits< ElementImp, ArrayImp > >
  {
  public:
    typedef ElementImp ElementType;

    typedef ArrayImp ArrayType;
    
    typedef ArrayDefaultTraits< ElementType, ArrayType > TraitsType;

  private:
    typedef ArrayDefault< ElementType, ArrayType > ThisType;
    typedef ArrayInterface< TraitsType > BaseType;

    using BaseType :: size;
    using BaseType :: asImp;

  public:
    typedef typename TraitsType :: IteratorType IteratorType;
    typedef typename TraitsType :: ConstIteratorType ConstIteratorType;

  public:
    inline ArrayImp& operator= ( const ElementType &element )
    {
      ArrayImp &imp = asImp();
      const unsigned int size = imp.size();
      for( unsigned int i = 0; i < size; ++i )
        imp[ i ] = element;
      return imp;
    }

    inline ConstIteratorType begin () const
    {
      return ConstIteratorType( asImp(), 0 );
    }

    inline IteratorType begin ()
    {
      return IteratorType( asImp(), 0 );
    }

    inline ConstIteratorType end () const
    {
      return ConstIteratorType( asImp(), size() );
    }

    inline IteratorType end ()
    {
      return IteratorType( asImp(), size() );
    }
  };



  template< class ElementImp >
  class ArrayWrapper
  : public ArrayDefault< ElementImp, ArrayWrapper< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef ArrayWrapper< ElementType > ThisType;
    typedef ArrayDefault< ElementType, ThisType > BaseType;

  private:
    const unsigned int size_;
    ElementType *elements_;

  public:
    inline ArrayWrapper ( unsigned int size, ElementType *elements )
    : size_( size ),
      elements_( elements )
    {
      assert( elements_ != NULL );
    }
    
    inline ThisType &operator= ( const ElementType &element )
    {
      return BaseType :: operator=( element );
    }

    inline const ElementType &operator[] ( unsigned int index ) const
    {
      assert( index < size_ );
      return elements_[ index ];
    }

    inline ElementType &operator[] ( unsigned int index )
    {
      assert( index < size_ );
      return elements_[ index ];
    }

    inline unsigned int size () const
    {
      return size_;
    }
  };



  template< class ElementImp, unsigned int arraysize >
  class StandardArray
  : public ArrayDefault< ElementImp, StandardArray< ElementImp, arraysize > >
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef StandardArray< ElementType, arraysize > ThisType;
    typedef ArrayDefault< ElementType, ThisType > BaseType;

  protected:
    ElementType elements_[ arraysize ];

  public:
    inline StandardArray ()
    {
    }

    inline StandardArray ( const ElementType &element )
    {
      for( unsigned int i = 0; i < arraysize; ++i )
        elements_[ i ] = element;
    }

    inline ThisType &operator= ( const ElementType &element )
    {
      return BaseType :: operator=( element );
    }
   
    inline const ElementType &operator[] ( unsigned int index ) const
    {
      assert( index < arraysize );
      return elements_[ index ];
    }

    inline ElementType &operator[] ( unsigned int index )
    {
      assert( index < arraysize );
      return elements_[ index ];
    }

    inline unsigned int size () const
    {
      return arraysize;
    }
  };



  template< class ElementImp, class ArrayAllocatorImp >
  class ArrayAllocatorInterface
  {
  public:
    typedef ElementImp ElementType;

    typedef ArrayAllocatorImp ArrayAllocatorType;

  private:
    typedef ArrayAllocatorInterface< ElementType, ArrayAllocatorType > ThisType;

  public:
    typedef ThisType ArrayAllocatorInterfaceType;

    typedef ElementType *ElementPtrType;

  public:
    inline ArrayAllocatorInterface ()
    {
      typedef CompileTimeChecker< Conversion< ArrayAllocatorType, ThisType > :: exists >
        __Array_Allocator_Implementation_Must_Be_Derived_From_Interface__;
    }
    
    inline void allocate ( unsigned int size,
                           ElementPtrType &array ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().allocate( size, array ) );
    }

    inline void free ( ElementPtrType &array ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().free( array ) );
    }

    inline void reallocate ( unsigned int oldSize,
                             unsigned int newSize,
                             ElementPtrType &array ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().reallocate( oldSize, newSize, array ) );
    }

  protected:
    inline const ArrayAllocatorType &asImp () const
    {
      return static_cast< const ArrayAllocatorType& >( *this );
    }

    inline ArrayAllocatorType &asImp ()
    {
      return static_cast< ArrayAllocatorType& >( *this );
    }
  };



  template< class ElementImp >
  class DefaultArrayAllocator
  : public ArrayAllocatorInterface< ElementImp, DefaultArrayAllocator< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef DefaultArrayAllocator< ElementType > ThisType;
    typedef ArrayAllocatorInterface< ElementType, ThisType > BaseType;

  public:
    typedef ElementType *ElementPtrType;
    
  public:
    inline void allocate ( unsigned int size,
                           ElementPtrType &array ) const
    {
      array = new ElementType[ size ];
      assert( array != 0 );
    }
  
    inline void free ( ElementPtrType &array ) const
    {
      delete[]( array );
      array = 0;
    }

    inline void reallocate ( unsigned int oldSize,
                             unsigned int newSize,
                             ElementPtrType &array ) const
    {
      ElementPtrType p = new ElementType[ newSize ];
      const unsigned int copySize = oldSize < newSize ? oldSize : newSize;
      for( unsigned int i = 0; i < copySize; ++i )
        p[ i ] = array[ i ];
      delete[]( array );
      array = p;
    }
  };



  template< class ElementImp,
            class ArrayAllocatorImp = DefaultArrayAllocator< ElementImp > >
  class DynamicArray
  : public ArrayDefault< ElementImp, DynamicArray< ElementImp, ArrayAllocatorImp > >
  {
  public:
    typedef ElementImp ElementType;

    typedef ArrayAllocatorImp ArrayAllocatorType;

  private:
    typedef DynamicArray< ElementType > ThisType;
    typedef ArrayDefault< ElementType, ThisType > BaseType;

  protected:
    ArrayAllocatorType allocator_;
    
    unsigned int size_;
    ElementType *elements_;

  public:
    inline DynamicArray ()
    {
      size_ = 0;
      allocator_.allocate( size_, elements_ );
    }
    
    inline DynamicArray ( unsigned int size )
    {
      size_ = size;
      allocator_.allocate( size_, elements_ );
    }
    
    inline DynamicArray ( unsigned int size, const ElementType &element )
    {
      size_ = size;
      allocator_.allocate( size_, elements_ );

      for( unsigned int i = 0; i < size_; ++i )
        elements_[ i ] = element;
    }


    inline ~DynamicArray ()
    {
      allocator_.free( elements_ );
    }

    inline ThisType& operator= ( const ElementType &element )
    {
      return BaseType :: operator=( element );
    }

    inline const ElementType& operator[] ( unsigned int index ) const
    {
      assert( index < size_ );
      return elements_[ index ];
    }

    inline ElementType& operator[] ( unsigned int index )
    {
      assert( index < size_ );
      return elements_[ index ];
    }

    inline void append ( const ElementType &element )
    {
      const unsigned int oldSize = size_;
      resize( oldSize + 1 );
      elements_[ oldSize ] = element;
    }

    template< class ArrayType >
    inline void append ( const ArrayInterface< typename ArrayType :: TraitsType > &array )
    {
      const unsigned int arraySize = array.size();
      
      const unsigned int oldSize = size_;
      resize( oldSize + arraySize );

      for( unsigned int i = 0; i < arraySize; ++i )
        elements_[ oldSize + i ] = array[ i ];
    }

    inline void resize ( unsigned int newSize )
    {
      allocator_.reallocate( size_, newSize, elements_ );
      size_ = newSize;
    }

    inline unsigned int size () const
    {
      return size_;
    }
  };

}

#endif
