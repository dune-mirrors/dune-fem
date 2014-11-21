#ifndef DUNE_FEM_METATWISTUTILITY_HH
#define DUNE_FEM_METATWISTUTILITY_HH

#include <cassert>

#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune
{

  namespace Fem
  {

    /** \brief MetaTwistUtility forwards the twist calls to the TwistUtility of the
     *         underlying HostTwistUtility.
     *
     *  \note The class Intersection implementation is assumed to have a method
     *        hostIntersection().
     */
    template< class HostTwistUtility >
    struct MetaTwistUtility
    {
      typedef HostTwistUtility  HostTwistUtilityType;
      typedef typename HostTwistUtilityType :: GridType  GridType;

      //! \brief return 0 for inner face
      template< class Intersection >
      static int twistInSelf ( const GridType & grid, const Intersection & intersection )
      {
        return HostTwistUtilityType::twistInSelf( grid, intersection.impl().hostIntersection() );
      }

      //! \brief return 0 for outer face
      template< class Intersection >
      static int twistInNeighbor ( const GridType & grid , const Intersection & intersection )
      {
        return HostTwistUtilityType::twistInNeighbor( grid, intersection.impl().hostIntersection() );
      }

      /** \brief return geometry type of inside or outside entity */
      template< class Intersection >
      static GeometryType elementGeometry ( const Intersection &intersection, const bool inside )
      {
        return HostTwistUtilityType::elementGeometry( intersection.impl().hostIntersection(), inside );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_METATWISTUTILITY_HH
