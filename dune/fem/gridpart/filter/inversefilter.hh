#ifndef DUNE_FEM_GRIDPART_FILTER_INVERSEFILTER_HH
#define DUNE_FEM_GRIDPART_FILTER_INVERSEFILTER_HH

namespace Dune
{

  namespace Fem
  {

    template< class FilterImp >
    class InverseFilter
    {
    public:
      typedef InverseFilter< FilterImp > ThisType;

      //! \brief entity types
      template< int cd >
      struct Codim
      {
        typedef typename FilterImp::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of entity with codim=0
      typedef typename Codim< 0 >::EntityType EntityType;

      //! \brief constructor
      InverseFilter ( const FilterImp & filter = FilterImp() )
      : filter_( filter )
      { }

      InverseFilter( const ThisType & ) = default;
      InverseFilter( ThisType && ) = default;

      ThisType &operator= ( const ThisType & ) = default;
      ThisType &operator= ( ThisType && ) = default;

      //! \brief returns true if the given entity of the pointer in the domain
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & entity ) const
      {
        return !filter().contains< cd >( entity );
      }

      //! \brief returns true if the given entity of the pointer in the domain
      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        return !filter().contains( entity );
      }

      //! \brief returns true if an intersection is interior
      //         (allows boundarys within a given domain)
      template< class Intersection >
      bool interiorIntersection ( const Intersection &intersection ) const
      {
        return contains( intersection.outside() );
      }

      //! \brief returns true if an intersection is a boundary intersection
      template< class Intersection >
      bool intersectionBoundary( const Intersection &intersection ) const
      {
        return filter().intersectionBoundary( intersection );
      }

      //! \brief returns the boundary id for an intersection
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection &intersection ) const
      {
        return filter().intersectionBoundaryId( intersection );
      }

      //! \brief returns true if for an intersection a neighbor exsits
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection &intersection ) const
      {
        return filter().intersectionNeighbor( intersection );
      }

    private:
      const FilterImp filter () const
      {
        return filter_;
      }

      FilterImp filter_;

    }; // end class InverseFilter

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTER_INVERSEFILTER_HH
