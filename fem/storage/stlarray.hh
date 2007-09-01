#ifndef DUNE_FEM_STLARRAY_HH
#define DUNE_FEM_STLARRAY_HH

#include <vector>

#include <dune/fem/storage/array.hh>

namespace Dune
{

  template< class ElementImp, class RealImp >
  class STLArrayIterator
  {
  public:
    typedef ElementImp ElementType;

  private:
    typedef STLArrayIterator< ElementType, RealImp > ThisType;

  protected:
    RealImp realIterator_;

  public:
    inline explicit STLArrayIterator ( RealImp realIterator )
    : realIterator_( realIterator )
    {
    }

    inline STLArrayIterator(  const ThisType &other )
    : realIterator_( other.realIterator )
    {
    }

    inline ThisType &operator= ( const ThisType &other )
    {
      realIterator_ = other.realIterator_;
    }

    inline ElementType &operator* ()
    {
      return *realIterator_;
    }

    inline ThisType &operator++ ()
    {
      ++realIterator_;
      return *this;
    }

    inline bool operator== ( const ThisType &other ) const
    {
      return realIterator_ == other.realIterator_;
    }

    inline bool operator!= ( const ThisType &other ) const
    {
      return realIterator_ != other.realIterator_;
    }
  };



  template< class ElementImp >
  class STLArray;



  template< class ElementImp >
  struct STLArrayTraits
  {
    typedef ElementImp ElementType;

    typedef STLArray< ElementType > ArrayType;

    typedef std :: vector< ElementType > stdVectorType;

    typedef STLArrayIterator< ElementType, typename stdVectorType :: iterator >
      IteratorType;
    typedef STLArrayIterator< const ElementType, typename stdVectorType :: const_iterator >
      ConstIteratorType;
  };



  template< class ElementImp >
  class STLArray
  : public ArrayInterface< STLArrayTraits< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;

    typedef STLArrayTraits< ElementType > Traits;

  private:
    typedef STLArray< ElementType > ThisType;
    typedef ArrayInterface< Traits > BaseType;

    typedef typename Traits :: stdVectorType stdVectorType;

  public:
    typedef typename Traits :: IteratorType IteratorType;
    typedef typename Traits :: ConstIteratorType ConstIteratorType;

  protected:
    stdVectorType vector_;

  public:
    inline explicit STLArray ( unsigned int size = 0 )
    : vector_( size )
    {
    }

    inline STLArray ( unsigned int size,
                      const ElementType &element )
    : vector_( size, element )
    {
    }

    inline STLArray ( const ThisType &other )
    : vector_( other.vector_ )
    {
    }

    inline const ElementType &operator[] ( unsigned int index ) const
    {
      return vector_[ index ];
    }

    inline ElementType &operator[] ( unsigned int index )
    {
      return vector_[ index ];
    }

    //! fill the array with copies of an element
    inline ThisType &assign ( const ElementType &element )
    {
      vector_.assign( size(), element );
      return *this;
    }

    //! copy another array to this one
    template< class T >
    inline ThisType &assign( const ArrayInterface< T > &other )
    {
      const unsigned int size = other.size();
      resize( size );
      for( unsigned int i = 0; i < size; ++i )
        vector_[ i ] = other[ i ];
      return *this;
    }

    inline void append ( const ElementType &element )
    {
      vector_.push_back( element );
    }

    template< class T >
    inline void append ( const ArrayInterface< T > &array )
    {
      const unsigned int arraySize = array.size();
      for( unsigned int i = 0; i < arraySize; ++i )
        append( array[ i ] );
    }

    inline void resize ( unsigned int newSize )
    {
      vector_.resize( newSize );
    }

    //! obtain begin iterator
    inline ConstIteratorType begin () const
    {
      return ConstIteratorType( vector_.begin() );
    }

    //! obtain begin iterator
    inline IteratorType begin ()
    {
      return IteratorType( vector_.begin() );
    }

    //! obtain end iterator
    inline ConstIteratorType end () const
    {
      return ConstIteratorType( vector_.end() );
    }

    //! obtain end iterator
    inline IteratorType end ()
    {
      return IteratorType( vector_.end() );
    }

    inline unsigned int size () const
    {
      return vector_.size();
    }
  };

}

#endif
