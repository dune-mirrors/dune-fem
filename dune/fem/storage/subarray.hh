#ifndef DUNE_FEM_SUBARRAY_HH
#define DUNE_FEM_SUBARRAY_HH

#include <type_traits>

#include <dune/fem/storage/vector.hh>

namespace Dune
{
  namespace Fem
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

      //! Maps an index onto another one
      unsigned int operator[] ( unsigned int index ) const
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
      static const bool v = std::is_convertible< IndexMapper, IndexMapperInterfaceType >::value;
    };



    // SubVector
    template< class BaseVectorImp, class IndexMapperImp >
    class SubVector
    : public VectorDefault< typename BaseVectorImp :: FieldType, SubVector< BaseVectorImp, IndexMapperImp > >
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

      BaseVectorType &baseVector_;
      const IndexMapperType &indexMapper_;

    public:
      SubVector( BaseVectorType &baseVector, const IndexMapperType &indexMapper )
      : baseVector_( baseVector ), indexMapper_( indexMapper )
      {
        assert( (unsigned int)baseVector_.size() == indexMapper_.range() );
      }

      SubVector ( const ThisType & ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;

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

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SUBARRAY_HH
