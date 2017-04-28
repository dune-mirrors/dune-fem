#ifndef DUNE_FEM_DOMAINFILTER_HH
#define DUNE_FEM_DOMAINFILTER_HH

#include <dune/fem/gridpart/filter/filter.hh>
#include <dune/fem/storage/dynamicarray.hh>

#include <dune/fem/misc/boundaryidprovider.hh>

namespace Dune
{
  namespace Fem
  {

    // forward declarations
    // --------------------

    template< class > class FilterDefaultImplementation;
    template< class , class > class DomainFilter;


    // DomainFilterTraits
    // ------------------------

    template< class GridPartImp, class DomainArrayImp >
    struct DomainFilterTraits
    {
      //! \brief grid part type
      typedef GridPartImp GridPartType;

      //! \brief array type
      typedef DomainArrayImp DomainArrayType;

      //! \brief filter type
      typedef DomainFilter< GridPartType, DomainArrayType > FilterType;

      //! \brief entity types
      template < int cd >
      struct Codim
      {
        //! \brief entity type for given codimension
        typedef typename GridPartType::template Codim< cd >::EntityType EntityType;
      };

      //! \brief entity type for codimension 0
      typedef typename Codim< 0 >::EntityType EntityType;
    };


    // DomainFilter
    // ------------------

    template< class GridPartImp, class DomainArrayImp = DynamicArray< int > >
    class DomainFilter
    : public FilterDefaultImplementation< DomainFilterTraits< GridPartImp, DomainArrayImp > >
    {
    public:
      //! \brief type of grid part
      typedef GridPartImp GridPartType;

      //! \brief type of array
      typedef DomainArrayImp DomainArrayType;

      //! \brief type of array field
      typedef typename DomainArrayType::value_type FieldType;

      //! \brief type of traits
      typedef DomainFilterTraits< GridPartType, DomainArrayType > Traits;

      //! \brief boundary id provider, specialized for each grid
      typedef BoundaryIdProvider< typename GridPartType::GridType >
        BoundaryIdProviderType;
    protected:
      //! \brief this type
      typedef DomainFilter< GridPartType, DomainArrayType > ThisType;

      //! \brief base type
      typedef FilterDefaultImplementation< Traits > BaseType;

    public:
      //! \brief type of the filter implementation
      typedef typename Traits::FilterType FilterType;

      //! \brief type of index set
      typedef typename GridPartType :: IndexSetType IndexSetType;

      template< int cd >
      struct Codim
      {
        typedef typename Traits::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of codim 0 entity
      typedef typename Traits::EntityType EntityType;

      //! \brief constructor
      DomainFilter ( const GridPartType & gridPart,
                     const DomainArrayType& tags,
                     const FieldType tag )
      : indexSet_( gridPart.indexSet() ),
        tags_( tags ),
        tag_( tag )
      {}

      DomainFilter ( const ThisType & ) = default;
      DomainFilter ( ThisType && ) = default;

      ThisType &operator= ( const ThisType & ) = default;
      ThisType &operator= ( ThisType && ) = default;

      //! \brief return false since all interior intersections should be skipped
      template< class Intersection >
      bool interiorIntersection( const Intersection & ) const
      {
        return false;
      }

      //! \brief returns true if the given entity has the correct tag
      //! for higher codims false is returned
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & entity ) const
      {
        if( cd == 0 )
        {
          return ( tag_ == tags_[ indexSet_.index( entity ) ] );
        }
        else
          return false;
      }

      //! \brief returns true if the given entity has the correct tag
      //! for higher codims false is returned
      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        return contains< Entity::codimension >( entity );
      }

      //! \brief returns true if an intersection is a boundary intersection
      template< class Intersection >
      bool intersectionBoundary( const Intersection & intersection ) const
      {
        return BoundaryIdProviderType::boundaryId( intersection );
      }

      //! \brief returns the boundary id for an intersection
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection & intersection ) const
      {
        return BoundaryIdProviderType::boundaryId( intersection );
      }

      //! \brief returns true if for an intersection a neighbor exsits
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection & intersection ) const
      {
        return intersection.neighbor();
      }

    protected:
      const IndexSetType& indexSet_;
      const DomainArrayType& tags_;
      const FieldType tag_;
    };

  }  // end namespace Fem

}  // end namespace Dune

#endif // #ifndef DUNE_FEM_DOMAINFILTER_HH
