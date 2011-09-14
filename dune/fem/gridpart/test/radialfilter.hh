#ifndef DUNE_FEM_GRIDPART_TEST_RADIALFILTER_HH
#define DUNE_FEM_GRIDPART_TEST_RADIALFILTER_HH

// dune-fem includes
#include <dune/fem/gridpart/filteredgridpart.hh>

namespace Dune
{
  namespace Fem
  {

    // forward declarations
    // --------------------

    template< class > struct DefaultFilterTraits;

    // RadialFilter
    // ------------

    /**! \brief example implementation; given center x and radius r, 
     *          filter is characteristic function of clos B_r( x )
     */
    template < class GridPartType >
    class RadialFilter
    {
    public:
      template< int cd >
      struct Codim
      {
        typedef typename GridPartType::template Codim< cd >::EntityType EntityType;
      };
      
      //! \brief type of entity
      typedef typename Codim< 0 >::EntityType EntityType;

      //! \brief geometry type
      typedef typename EntityType::Geometry GeometryType;

      //! \brief coordinate type
      typedef typename GeometryType::GlobalCoordinate GlobalCoordinateType;

      // traits
      typedef DefaultFilterTraits< RadialFilter< GridPartType > > Traits;

      //! \brief constructor
      RadialFilter( const GlobalCoordinateType & center, 
                    const double radius ) 
      : center_( center ),
        radius_( radius )
      { }

      //! \brief check whether barycenter is inside of circle 
      template< int cc >
      bool contains ( const EntityType & e ) const 
      {
        if( cc != 0 )
          DUNE_THROW( InvalidStateException, "RadialFilter::contains only available for codim 0 entities" );
        const GeometryType & geometry = e.geometry();
        double dist = (geometry.center() - center_).two_norm();
        return (dist <= radius_);
      }

      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        static const int cc = Entity::codimension;
        return contains< cc >( entity );
      }
      
      //! \brief return what boundary id we have in case of boundary intersection 
      //         which is either it.boundary == true or contains (it.ouside()) == false 
      //         so here true is a good choice 
      template < class IntersectionIteratorType >
      inline bool intersectionBoundary( const IntersectionIteratorType & it ) const 
      {
        return true;
      }
      //! \brief return what boundary id we have in case of boundary intersection 
      //          which is either it.boundary == true or contains (it.ouside()) == false 
      template < class IntersectionIteratorType >
      inline int intersectionBoundaryId(const IntersectionIteratorType & it) const
      {
        return 1;
      }

      //! \brief if contains() is true then we have an interior entity 
      template <class IntersectionIteratorType>
      inline bool intersectionNeighbor( const IntersectionIteratorType & it ) const
      {
        return true;
      }

    private:
      const GlobalCoordinateType center_;
      const double radius_;

    }; // end RadialFilter

  }  // end namespace Fem

}  // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_RADIALFILTER_HH

