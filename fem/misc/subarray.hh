#ifndef DUNE_FEM_SUBARRAY_HH
#define DUNE_FEM_SUBARRAY_HH

#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/misc/array.hh>

namespace Dune
{

  //! Abstract index mapper interface
  template< class Imp >
  class IndexMapperInterface
  {
  private:
    typedef IndexMapperInterface< Imp > ThisType;

  public:
    typedef ThisType IndexMapperInterfaceType;
      
  public:
    //! Maps an index onto another one
    inline const unsigned int &operator[]( unsigned int index ) const
    {
      return asImp().operator[]( index );
    }
    
    //! Returns the map's range
    inline unsigned int range () const
    {
      return asImp().size();
    }

    //! Returns the map's size
    inline unsigned int size () const
    {
      return asImp().size();
    }

  private:
    //! Barton-Nackman trick
    inline const Imp &asImp () const
    {
      return static_cast< const Imp& >( *this );
    }

    //! Barton-Nackman trick
    inline Imp &asImp ()
    {
      return static_cast< Imp& >( *this );
    }
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
      typedef ArrayInterface< typename BaseArrayType :: TraitsType >
        BaseArrayInterfaceType;
      typedef CompileTimeChecker
        < Conversion< BaseArrayType, BaseArrayInterfaceType > :: exists >
        __BaseArrayType_Must_Be_Derived_From_ArrayInterface__;

      typedef IndexMapperInterface< IndexMapperType >
        IndexMapperInterfaceType;
      typedef CompileTimeChecker
        < Conversion< IndexMapperType, IndexMapperInterfaceType > :: exists >
        __IndexMapperType_Must_Be_Derived_From_IndexMapperInterface__;        

      assert( baseArray_.size() == indexMapper_.range() );
    }

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

}

#endif
