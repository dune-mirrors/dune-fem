#ifndef DUNE_FEM_ARRAY_HH
#define DUNE_FEM_ARRAY_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  template< class ElementImp, class ArrayImp >
  class ArrayInterface
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef ArrayInterface< ElementType, ArrayImp > ThisType;

  public:
    inline ArrayImp& operator= ( const ElementType &element )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator=( element ) );
      return asImp();
    }
    
    inline const ElementType& operator[] ( unsigned int index ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
    }

    inline ElementType& operator[] ( unsigned int index )
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
      return asImp()[ index ];
    }

    inline unsigned int size () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
      return asImp().size();
    }
    
  protected:
    inline const ArrayImp& asImp () const
    {
      return static_cast< const ArrayImp& >( *this );
    }

    inline ArrayImp& asImp ()
    {
      return static_cast< ArrayImp& >( *this );
    }
  };


  
  template< class ElementImp, class ArrayImp >
  class ArrayDefault
  : public ArrayInterface< ElementImp, ArrayImp >
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef ArrayDefault< ElementType, ArrayImp > ThisType;
    typedef ArrayInterface< ElementType, ArrayImp > BaseType;

  public:
    inline ArrayImp& operator= ( const ElementType &element )
    {
      ArrayImp &imp = this->asImp();
      const unsigned int size = imp.size();
      for( unsigned int i = 0; i < size; ++i )
        imp[ i ] = element;
      return imp;
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

    inline unsigned int size () const
    {
      return size_;
    }
  };



  template< class ElementImp >
  class DynamicArray
  : public ArrayDefault< ElementImp, DynamicArray< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef DynamicArray< ElementType > ThisType;
    typedef ArrayDefault< ElementType, ThisType > BaseType;

  private:
    static unsigned int defaultSize_;
    
    const unsigned int size_;
    ElementType *elements_;

  public:
    inline DynamicArray ()
    : size_( defaultSize_ ),
      elements_( new ElementType[ size_ ] )
    {
      assert( elements_ != NULL );
    }
    
    inline DynamicArray ( unsigned int size )
    : size_( size ),
      elements_( new ElementType[ size_ ] )
    {
      assert( elements_ != NULL );
    }
    
    inline DynamicArray ( unsigned int size, const ElementType &element )
    : size_( size ),
      elements_( new ElementType[ size_ ] )
    {
      assert( elements_ != NULL );
      for( unsigned int i = 0; i < size_; ++i )
        elements_[ i ] = element;
    }


    inline ~DynamicArray ()
    {
      delete[] elements_;
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

    inline static unsigned int defaultSize ( unsigned int size )
    {
      unsigned int current = defaultSize_;
      defaultSize_ = size;
      return current;
    }

    inline unsigned int size () const
    {
      return size_;
    }
  };
  
  template< class ElementType >
  unsigned int DynamicArray< ElementType > :: defaultSize_ = 0;
}

#endif
