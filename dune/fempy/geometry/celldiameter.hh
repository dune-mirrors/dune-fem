#ifndef DUNE_FEMPY_GEOMETRY_CELLDIAMETER_HH
#define DUNE_FEMPY_GEOMETRY_CELLDIAMETER_HH

#include <cassert>
#include <dune/common/math.hh>
#include <dune/fempy/geometry/edgelength.hh>

#if __GNUC__ >= 13
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-reference"
#endif

namespace Dune
{

  namespace FemPy
  {
    // cellDiameter
    // ------------
    //
    // computes diam(E) = sup_{x,y \in E} || x - y ||
    //
    template< class Geometry >
    inline static typename Geometry::ctype
    cellDiameter ( const Geometry &geometry )
    {
      typedef typename Geometry::ctype ctype;
      if constexpr ( Geometry::mydimension > 0)
      {
        // for simplices the maxEdgeLength is the same as cellDiameter
        const auto type = geometry.type();
        if( type.isSimplex() )
        {
          return maxEdgeLength( geometry );
        }
        else if ( type.isCube() ) // 2d quad and 3d hexa
        {
          assert( type.isQuadrilateral() || type.isHexahedron() );
          ctype diam = 0.;
          const int n = Dune::power( 2, Geometry::mydimension-1 );
          const int last = geometry.corners() - 1;
          for( int i=0; i<n; ++i )
          {
            auto dist = geometry.corner(i);
            dist -= geometry.corner(last - i);
            diam = std::max( diam, dist.two_norm2() );
          }
          return std::sqrt( diam );
        }
        else // all other types
        {
          ctype diam = 0.;
          const int corners = geometry.corners();
          for( int i=0; i<corners; ++i)
          {
            const auto x = geometry.corner( i );
            for( int j=i+1; j<corners; ++j )
            {
              auto dist = geometry.corner( j );
              dist -= x;
              diam = std::max( diam, dist.two_norm2() );
            }
          }
          return std::sqrt( diam );
        }
      }
      return 1;
    }

  } // namespace FemPy

} // namespace Dune

#if __GNUC__ >= 13
#pragma GCC diagnostic pop
#endif

#endif // #ifndef DUNE_FEMPY_GEOMETRY_EDGELENGTH_HH
