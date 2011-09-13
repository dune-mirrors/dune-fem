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

    template < class, class > struct DefaultFilterTraits;
    template< class > class FilterDefaultImplementation;

    // RadialFilter
    // ------------

    /**! \brief example implementation; given center x and radius r, 
     *          filter is characteristic function of clos B_r( x )
     */
    template < class GridPartType >
    class RadialFilter
    : public FilterDefaultImplementation< DefaultFilterTraits< RadialFilter< GridPartType >, GridPartType > >
    {
      // type of this
      typedef RadialFilter< GridPartType > ThisType;

      // traits
      typedef DefaultFilterTraits< RadialFilter< GridPartType >, GridPartType > Traits;

      // base type
      typedef FilterDefaultImplementation< Traits > BaseType;

      // geometry type
      typedef typename BaseType::EntityType::Geometry GeometryType;

    public:
      //! \brief type of entity
      typedef typename BaseType::EntityType EntityType;

      //! \brief type of entity pointer
      typedef typename BaseType::EntityPointerType EntityPointerType;

      //! \brief coordinate type
      typedef typename GeometryType::GlobalCoordinate GlobalCoordinateType;

      //! \brief constructor
      RadialFilter( const GlobalCoordinateType & center, 
                    const double radius ) 
      : center_( center ),
        radius_( radius )
      { }

      //! \brief check whether barycenter is inside of circle 
      inline bool contains ( const EntityPointerType & ep ) const
      {
        return contains( *ep );
      }
      
      //! \brief check whether barycenter is inside of circle 
      inline bool contains ( const EntityType & e ) const 
      {
        const GeometryType & geometry = e.geometry();
        double dist = (geometry.center() - center_).two_norm();
        return (dist <= radius_);
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

      //! \brief create default radial filter 
      static ThisType createObject( const GridPartType & gridPart )
      {
        std::cerr << "Warning, creating default readial filter! " << std::endl;
        GlobalCoordinateType  center( 0 );
        double radius = 0.5;
        return ThisType( center, radius );
      }

    private:
      const GlobalCoordinateType center_;
      const double radius_;

    }; // end RadialFilter

  }  // end namespace Fem

}  // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_RADIALFILTER_HH

