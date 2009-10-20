#ifndef DUNE_FEM_SUBARRAY_HH
#define DUNE_FEM_SUBARRAY_HH

#include <dune/common/typetraits.hh>

#include <dune/fem/storage/array.hh>
#include <dune/fem/storage/vector.hh>

namespace Dune
{

  //! Abstract index mapper interface
  template< class IM >
  class IndexMapperInterface
  : public BartonNackmanInterface< IndexMapperInterface< IM >, IM >
  {
    typedef IndexMapperInterface< IM > ThisType;
    typedef BartonNackmanInterface< ThisType, IM > BaseType;

  public:
    //! type of the implementation (Barton-Nackman)
    typedef IM IndexMapperType;
    
    //! type of the interface
    typedef ThisType IndexMapperInterfaceType;

  public:
    //! Maps an index onto another one
    const unsigned int operator[] ( unsigned int index ) const
    {
      return asImp().operator[]( index );
    }
    
    //! Returns the map's range
    unsigned int range () const
    {
      return asImp().range();
    }

    //! Returns the map's size
    unsigned int size () const
    {
      return asImp().size();
    }

  protected:
    using BaseType::asImp;
  };
  
  
  
  template< class IndexMapper >
  struct SupportsIndexMapperInterface
  {
    typedef IndexMapperInterface< IndexMapper > IndexMapperInterfaceType;
    static const bool v = Conversion< IndexMapper, IndexMapperInterfaceType >::exists;
  };



  // SubArray
  template< class BaseArrayImp, class IndexMapperImp >
  class SubArray
  : public ArrayDefault< typename BaseArrayImp :: ElementType,
                         SubArray< BaseArrayImp, IndexMapperImp > >
  {
  public:
    //! type of the base array
    typedef BaseArrayImp BaseArrayType;

    //! type of the index mapper
    typedef IndexMapperImp IndexMapperType;

    //! type of array elements
    typedef typename BaseArrayType :: ElementType ElementType;
      
  private:
    typedef SubArray< BaseArrayType, IndexMapperType > ThisType;
    typedef ArrayDefault< ElementType, ThisType > BaseType;
      
  private:
    BaseArrayType &baseArray_;
    const IndexMapperType &indexMapper_;
      
  public:
    SubArray( BaseArrayType &baseArray, const IndexMapperType &indexMapper )
    : baseArray_( baseArray ),
      indexMapper_( indexMapper )
    {
      dune_static_assert( SupportsArrayInterface< BaseArrayType >::v, "SubArray can only wrap arrays." );
      dune_static_assert( SupportsIndexMapperInterface< IndexMapperType >::v, "Invalid index mapper." );
      assert( baseArray_.size() == indexMapper_.range() );
    }

    SubArray ( const ThisType &other )
    : baseArray_( other.baseArray_ ),
      indexMapper_( other.indexMapper_ )
    {}

  private:
    ThisType &operator= ( const ThisType &other );

  public:
    const ElementType &operator[] ( unsigned int index ) const
    {
      return baseArray_[ indexMapper_[ index ] ];
    }

    ElementType &operator[] ( unsigned int index )
    {
      return baseArray_[ indexMapper_[ index ] ];
    }

    unsigned int size () const
    {
      return indexMapper_.size();
    }
  };



  // SubVector
  template< class BaseVectorImp, class IndexMapperImp >
  class SubVector
  : public VectorDefault< typename BaseVectorImp :: FieldType,
                          SubVector< BaseVectorImp, IndexMapperImp > >
  {
  public:
    //! type of the base array
    typedef BaseVectorImp BaseVectorType;

    //! type of the index mapper
    typedef IndexMapperImp IndexMapperType;

    //! type of array elements
    typedef typename BaseVectorType :: FieldType FieldType;
      
  private:
    typedef SubVector< BaseVectorType, IndexMapperType > ThisType;
    typedef VectorDefault< FieldType, ThisType > BaseType;
      
  private:
    BaseVectorType &baseVector_;
    const IndexMapperType &indexMapper_;

  public:
    SubVector( BaseVectorType &baseVector, const IndexMapperType &indexMapper )
    : baseVector_( baseVector ),
      indexMapper_( indexMapper )
    {
      // used for StaticArray which only implements parts of the VectorInterface (axpy <-> addScaled for example)
      // but SubVector obly requires parts of the VectorInterface (operator[]...)
      // dune_static_assert( SupportsVectorInterface< BaseVectorType >::v, "SubVector can only wrap vectors." );
      // dune_static_assert( SupportsIndexMapperInterface< IndexMapperType >::v, "Invalid index mapper." );

      assert( (unsigned int)baseVector_.size() == indexMapper_.range() );
    }
    

  private:
    SubVector ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    const FieldType &operator[] ( unsigned int index ) const
    {
      return baseVector_[ indexMapper_[ index ] ];
    }

    FieldType &operator[] ( unsigned int index )
    {
      return baseVector_[ indexMapper_[ index ] ];
    }

    unsigned int size () const
    {
      return indexMapper_.size();
    }
  };
  
}

#endif // #ifndef DUNE_FEM_SUBARRAY_HH
