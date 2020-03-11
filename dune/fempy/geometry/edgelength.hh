#ifndef DUNE_FEMPY_GEOMETRY_EDGELENGTH_HH
#define DUNE_FEMPY_GEOMETRY_EDGELENGTH_HH

#include <limits>

#include <dune/geometry/referenceelements.hh>

namespace Dune
{

  namespace FemPy
  {

    // edgeLength
    // ----------

    template< class Geometry, class RefEmbedding >
    inline static typename Geometry::ctype edgeLength ( const Geometry &geometry, const RefEmbedding &refEmbedding )
    {
      typedef typename Geometry::ctype ctype;

      // refEmbedding is affine, so the tangent vector is constant
      const auto &refTangent = refEmbedding.jacobianTransposed( { ctype( 1 ) / ctype( 2 ) } )[ 0 ];

      // Use one-point quadrature to compute edge length
      typename Geometry::GlobalCoordinate tangent( 0 );
      geometry.jacobianTransposed( refEmbedding.center() ).umtv( refTangent, tangent );
      return tangent.two_norm();
    }

    template< class Geometry >
    inline static typename Geometry::ctype edgeLength ( const Geometry &geometry, int i )
    {
      const int dim = Geometry::mydimension;
      const auto &refElement = ReferenceElements< typename Geometry::ctype, dim >::general( geometry.type() );
      return edgeLength( geometry, refElement.template geometry< dim-1 >( i ) );
    }



    // maxEdgeLength
    // -------------

    template< class Geometry >
    inline static typename Geometry::ctype
    maxEdgeLength ( const Geometry &geometry )
    {
      if constexpr ( Geometry::mydimension > 0)
      {
        typedef typename Geometry::ctype ctype;
        const int dim = Geometry::mydimension;

        const auto &refElement = ReferenceElements< ctype, dim >::general( geometry.type() );

        ctype maxEdgeLength( 0 );
        for( int i = 0, n = refElement.size( dim-1 ); i < n; ++i )
          maxEdgeLength = std::max( maxEdgeLength, edgeLength( geometry, refElement.template geometry< dim-1 >( i ) ) );

        return maxEdgeLength;
      }
      return 1;
    }


    // minEdgeLength
    // -------------

    template< class Geometry >
    inline static typename Geometry::ctype minEdgeLength ( const Geometry &geometry )
    {
      typedef typename Geometry::ctype ctype;
      const int dim = Geometry::mydimension;

      const auto &refElement = ReferenceElements< ctype, dim >::general( geometry.type() );

      ctype minEdgeLength = std::numeric_limits< ctype >::max();
      for( int i = 0, n = refElement.size( dim-1 ); i < n; ++i )
        minEdgeLength = std::min( minEdgeLength, edgeLength( geometry, refElement.template geometry< dim-1 >( i ) ) );

      return minEdgeLength;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GEOMETRY_EDGELENGTH_HH
