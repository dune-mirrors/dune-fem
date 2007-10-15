#ifndef DUNE_FEM_ARRAY_HH
#define DUNE_FEM_ARRAY_HH

#include <cassert>

#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/storage/arrayallocator.hh>

namespace Dune
{

  /** \class ArrayInterface
   *  \brief abstract array interface
   */
  template< class TraitsImp >
  class ArrayInterface
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

  private:
    typedef ArrayInterface< Traits > ThisType;

  public:
    //! type of this interface
    typedef ThisType ArrayInterfaceType;
    
    //! type of the implementation (Barton-Nackman) 
    typedef typename Traits :: ArrayType ArrayType;

    //! type of the array elements
    typedef typename Traits :: ElementType ElementType;

    //! type of constant iterator
    typedef typename Traits :: ConstIteratorType ConstIteratorType;

    //! type of (non-constant) iterator
    typedef typename Traits :: IteratorType IteratorType;

  public:
    inline ArrayInterface ()
    {
      // make sure the implementation is derived from this interface
      typedef CompileTimeChecker< Conversion< ArrayType, ThisType > :: exists >
        __Array_Implementation_Must_Be_Derived_From_Interface__;
    }

    /** \brief access an array element
     *
     *  \param[in]  index  index of the array element to access
     *  
     *  \returns a const reference to the array element
     */
    inline const ElementType &operator[] ( unsigned int index ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
    }
    
    /** \brief access an array element
     *
     *  \param[in]  index  index of the array element to access
     *  
     *  \returns a reference to the array element
     */
    inline ElementType &operator[] ( unsigned int index )
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
    }
   
    /** \brief fill the array with copies of an element
     *
     *  \param[in]  element  element wich shall be copied into every array
     *                       entry
     *
     *  \returns a reference to this array
     */
    inline ArrayType &assign ( const ElementType &element )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( element ) );
      return asImp();
    }

    /** \brief copy another array to this one
     *
     *  Copies the data from another array to this one. Both arrays must be of
     *  the same size.
     *
     *  \param[in]  other  array to copy
     *
     *  \returns a reference this this array
     */
    template< class T >
    inline ArrayType &assign( const ArrayInterface< T > &other )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( other ) );
      return asImp();
    }
 
    /** \brief obtain begin iterator
     *
     *  \returns an iterator pointing to the first array element
     */
    inline ConstIteratorType begin () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
      return asImp().begin();
    }

    /** \brief obtain begin iterator
     *
     *  \returns an iterator pointing to the first array element
     */
    inline IteratorType begin ()
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
      return asImp().begin();
    }

    /** \brief obtain end iterator
     *
     *  \returns an iterator pointing behind the last array element
     */
    inline ConstIteratorType end () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
      return asImp().end();
    }

    /** \brief obtain end iterator
     *
     *  \returns an iterator pointing behind the last array element
     */
    inline IteratorType end ()
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
      return asImp().end();
    }

    /** obtain the size of the array
     *
     *  \returns the size of the array
     */
    inline unsigned int size () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
      return asImp().size();
    }
    
  protected:
    // Barton-Nackman trick
    inline const ArrayType &asImp () const
    {
      return static_cast< const ArrayType& >( *this );
    }

    // Barton-Nackman trick
    inline ArrayType &asImp ()
    {
      return static_cast< ArrayType& >( *this );
    }
  };



  // helper structure to make sure an array is derived from the interface
  template< class ArrayType >
  struct CheckArrayInterface
  {
    typedef ArrayInterface< typename ArrayType :: Traits > ArrayInterfaceType;

    typedef CompileTimeChecker< Conversion< ArrayType, ArrayInterfaceType > :: exists >
      CheckerType;
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
  struct ArrayDefaultTraits
  {
    typedef ElementImp ElementType;

    typedef ArrayImp ArrayType;

    typedef ArrayDefaultIterator< ElementType, ArrayType > IteratorType;
    
    typedef ArrayDefaultIterator< const ElementType, const ArrayType > ConstIteratorType;
  };


  
  /** \class ArrayDefault
   *  \brief default implementation of the ArrayInterface
   */
  template< class ElementImp, class ArrayImp >
  class ArrayDefault
  : public ArrayInterface< ArrayDefaultTraits< ElementImp, ArrayImp > >
  {
  public:
    //! type of the array elements
    typedef ElementImp ElementType;

    //! type of the implementation (Barton-Nackman) 
    typedef ArrayImp ArrayType;
    
    //! type of the traits
    typedef ArrayDefaultTraits< ElementType, ArrayType > Traits;

  private:
    typedef ArrayDefault< ElementType, ArrayType > ThisType;
    typedef ArrayInterface< Traits > BaseType;

  public:
    using BaseType :: size;

  protected:
    using BaseType :: asImp;

  public:
    //! type of constant iterator
    typedef typename Traits :: ConstIteratorType ConstIteratorType;

    //! type of (non-constant) iterator
    typedef typename Traits :: IteratorType IteratorType;

  public:
    /** \copydoc Dune::ArrayInterface::assign(const ElementType &element) */
    inline ArrayType &assign ( const ElementType &element )
    {
      ArrayType &imp = asImp();
      const unsigned int size = imp.size();
      for( unsigned int i = 0; i < size; ++i )
        imp[ i ] = element;
      return imp;
    }
    
    /** \copydoc Dune::ArrayInterface::assign(const ArrayInterface<T> &other) */
    template< class T >
    inline ArrayType &assign( const ArrayInterface< T > &other )
    {
      ArrayType &imp = asImp();
      const unsigned int size = imp.size();
      assert( size == other.size() );
      for( unsigned int i = 0; i < size; ++i )
        imp[ i ] = other[ i ];
      return imp;
    }

    /** \copydoc Dune::ArrayInterface::begin() const */
    inline ConstIteratorType begin () const
    {
      return ConstIteratorType( asImp(), 0 );
    }

    /** \copydoc Dune::ArrayInterface::begin() */
    inline IteratorType begin ()
    {
      return IteratorType( asImp(), 0 );
    }

    /** \copydoc Dune::ArrayInterface::end() const */
    inline ConstIteratorType end () const
    {
      return ConstIteratorType( asImp(), size() );
    }

    /** \copydoc Dune::ArrayInterface::end() */
    inline IteratorType end ()
    {
      return IteratorType( asImp(), size() );
    }
  };



  /** \class ArrayWrapper
   *  \brief implementation of the ArrayInterface wrapping a pointer to an
   *         array of elements
   */
  template< class ElementImp >
  class ArrayWrapper
  : public ArrayDefault< ElementImp, ArrayWrapper< ElementImp > >
  {
  public:
    //! type of the array elements
    typedef ElementImp ElementType;

  private:
    typedef ArrayWrapper< ElementType > ThisType;
    typedef ArrayDefault< ElementType, ThisType > BaseType;

  private:
    const unsigned int size_;
    ElementType *elements_;

  public:
    /** \brief create an ArrayWrapper from a size and a pointer */
    inline ArrayWrapper ( unsigned int size, ElementType *elements )
    : size_( size ),
      elements_( elements )
    {
      assert( elements_ != NULL );
    }
   
    /** \copydoc Dune::ArrayInterface::operator[](unsigned int index) const */
    inline const ElementType &operator[] ( unsigned int index ) const
    {
      assert( index < size_ );
      return elements_[ index ];
    }

    /** \copydoc Dune::ArrayInterface::operator[](unsigned int index) */
    inline ElementType &operator[] ( unsigned int index )
    {
      assert( index < size_ );
      return elements_[ index ];
    }

    /** \copydoc Dune::ArrayInterface::size */
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

    inline explicit StandardArray ( const ElementType &element )
    {
      assign( element );
    }

    inline StandardArray ( const ThisType &other )
    {
      assign( other );
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



  template< class ElementImp,
            template< class > class ArrayAllocatorImp = DefaultArrayAllocator >
  class DynamicArray
  : public ArrayDefault< ElementImp, DynamicArray< ElementImp, ArrayAllocatorImp > >
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef DynamicArray< ElementType, ArrayAllocatorImp > ThisType;
    typedef ArrayDefault< ElementType, ThisType > BaseType;

  protected:
    typedef ArrayAllocatorImp< ElementType > ArrayAllocatorType;
    
    typedef typename ArrayAllocatorType :: ElementPtrType ElementPtrType;

  protected:
    ArrayAllocatorType allocator_;
    
    unsigned int size_;
    ElementPtrType elements_;

  public:
    inline explicit DynamicArray ( unsigned int size = 0 )
    : allocator_()
    {
      size_ = size;
      allocator_.allocate( size_, elements_ );
    }

    inline explicit DynamicArray ( const ArrayAllocatorType &arrayAllocator,
                                   unsigned int size = 0 )
    : allocator_( arrayAllocator )
    {
      size_ = size;
      allocator_.allocate( size_, elements_ );
    }
    
    inline DynamicArray ( unsigned int size,
                          const ElementType &element )
    : allocator_()
    {
      size_ = size;
      allocator_.allocate( size_, elements_ );
      assign( element );
    }

    inline DynamicArray ( const ArrayAllocatorType &arrayAllocator,
                          unsigned int size,
                          const ElementType defaultElement )
    : allocator_( arrayAllocator )
    {
      size_ = size;
      allocator_.allocate( size_, elements_ );
      assign( defaultElement );
    }
    
    inline DynamicArray ( const ThisType &other )
    : allocator_( other.allocator_ )
    {
      size_ = other.size_;
      allocator_.allocate( size_, elements_ );
      assign( other );
    }

    inline ~DynamicArray ()
    {
      allocator_.free( elements_ );
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

    inline void append ( const ElementType element )
    {
      const unsigned int oldSize = size_;
      resize( oldSize + 1 );
      elements_[ oldSize ] = element;
    }

    template< class T >
    inline void append ( const ArrayInterface< T > &array )
    {
      const unsigned int arraySize = array.size();
      
      const unsigned int oldSize = size_;
      resize( oldSize + arraySize );

      for( unsigned int i = 0; i < arraySize; ++i )
        elements_[ oldSize + i ] = array[ i ];
    }

    inline void resize ( unsigned int newSize )
    {
      const unsigned int oldSize = size_;
      if( newSize == oldSize )
        return;

      allocator_.reallocate( oldSize, newSize, elements_ );
      size_ = newSize;
    }

    inline void resize ( unsigned int newSize,
                         const ElementType defaultElement )
    {
      const unsigned int oldSize = size_;
      if( newSize == oldSize )
        return;

      allocator_.reallocate( oldSize, newSize, elements_ );
      size_ = newSize;
      for( unsigned int i = oldSize; i < newSize; ++i )
        elements_[ i ] = defaultElement;
    }

    inline unsigned int size () const
    {
      return size_;
    }
  };

}

#include "array_inline.hh"

#endif
