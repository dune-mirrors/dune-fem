#ifndef DUNE_FEM_SUBARRAY_HH
#define DUNE_FEM_SUBARRAY_HH

#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/storage/array.hh>
#include <dune/fem/storage/vector.hh>

namespace Dune
{

  //! Abstract index mapper interface
  template< class IndexMapperImp >
  class IndexMapperInterface
  {
  public:
    //! type of the implementation (Barton-Nackman)
    typedef IndexMapperImp IndexMapperType;
    
  private:
    typedef IndexMapperInterface< IndexMapperType > ThisType;

  public:
    //! type of the interface
    typedef ThisType IndexMapperInterfaceType;

  public:
    //! Maps an index onto another one
    inline const unsigned int operator[] ( unsigned int index ) const
    {
      return asImp().operator[]( index );
    }
    
    //! Returns the map's range
    inline unsigned int range () const
    {
      return asImp().range();
    }

    //! Returns the map's size
    inline unsigned int size () const
    {
      return asImp().size();
    }

  private:
    //! Barton-Nackman trick
    inline const IndexMapperType &asImp () const
    {
      return static_cast< const IndexMapperType& >( *this );
    }

    //! Barton-Nackman trick
    inline IndexMapperType &asImp ()
    {
      return static_cast< IndexMapperType& >( *this );
    }
  };
  
  
  
  template< class IndexMapperType >
  struct CheckIndexMapperInterface
  {
    typedef IndexMapperInterface< IndexMapperType > IndexMapperInterfaceType;

    typedef CompileTimeChecker< Conversion< IndexMapperType, IndexMapperInterfaceType > :: exists >
      CheckerType;
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
    inline SubArray( BaseArrayType &baseArray,
                     const IndexMapperType &indexMapper )
    : baseArray_( baseArray ),
      indexMapper_( indexMapper )
    {
      typedef CheckArrayInterface< BaseArrayType > __CheckBaseArrayType__;
      typedef CheckIndexMapperInterface< IndexMapperType > __CheckIndexMapperType__;
      assert( baseArray_.size() == indexMapper_.range() );
    }

    inline SubArray ( const ThisType &other )
    : baseArray_( other.baseArray_ ),
      indexMapper_( other.indexMapper_ )
    {}

  private:
    ThisType &operator= ( const ThisType &other );

  public:
    inline const ElementType &operator[] ( unsigned int index ) const
    {
      return baseArray_[ indexMapper_[ index ] ];
    }

    inline ElementType &operator[] ( unsigned int index )
    {
      return baseArray_[ indexMapper_[ index ] ];
    }

    inline unsigned int size() const
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
    inline SubVector( BaseVectorType &baseVector,
                      const IndexMapperType &indexMapper )
    : baseVector_( baseVector ),
      indexMapper_( indexMapper )
    {
      typedef CheckVectorInterface< BaseVectorType > __CheckBaseVectorType__;
      typedef CheckIndexMapperInterface< IndexMapperType > __CheckIndexMapperType__;

      assert( (unsigned int)baseVector_.size() == indexMapper_.range() );
    }
    

  private:
    inline SubVector ( const ThisType &other );
    // : baseVector_( other.baseVector_ ),
    //  indexMapper_( other.indexMapper_ )
    //{}
    ThisType &operator= ( const ThisType &other );

  public:
    inline const FieldType &operator[] ( unsigned int index ) const
    {
      return baseVector_[ indexMapper_[ index ] ];
    }

    inline FieldType &operator[] ( unsigned int index )
    {
      return baseVector_[ indexMapper_[ index ] ];
    }

    inline unsigned int size() const
    {
      return indexMapper_.size();
    }
  };
  
}

#endif
