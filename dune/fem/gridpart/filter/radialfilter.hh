#ifndef DUNE_FEM_GRIDPART_FILTER_RADIALFILTER_HH
#define DUNE_FEM_GRIDPART_FILTER_RADIALFILTER_HH

// dune-common includes
#include <dune/common/fvector.hh>

namespace Dune
{

  namespace Fem
  {

    // RadialFilter
    // ------------

    /**! \brief example implementation; given center x and radius r,
     *          filter is characteristic function of clos B_r( x )
     */
    template < typename ct, int dimw >
    class RadialFilter
    {
    public:
      //! \brief export template parameter
      typedef ct ctype;

      //! \brief export dimension
      static const int dimensionworld = dimw;

      //! \brief coordinate type
      typedef Dune::FieldVector< ct, dimw > GlobalCoordinateType;

      //! \brief constructor
      RadialFilter( const GlobalCoordinateType & center,
                    const ctype radius )
      : center_( center ),
        radius_( radius )
      { }

      RadialFilter ()
      : center_( 0 ),
        radius_( 0.25 )
      { }

      //! \brief check whether entity center is inside of circle
      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        static const int cc = Entity::codimension;
        if( cc != 0 )
          DUNE_THROW( InvalidStateException, "RadialFilter::contains only available for codim 0 entities" );
        ctype dist = (entity.geometry().center() - center_).two_norm();
        return (dist > radius_);
      }

      //! \brief default implementation returns contains from neighbor
      template< class Intersection >
      bool interiorIntersection( const Intersection &intersection ) const
      {
        typedef typename Intersection::Entity EntityType;
        const EntityType& outside(intersection.outside());
        return contains( outside );
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
        return false;
      }

    private:
      const GlobalCoordinateType center_;
      const ctype radius_;

    }; // end RadialFilter

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTER_RADIALFILTER_HH

