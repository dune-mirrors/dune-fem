#include <config.h>

#include <dune/grid/common/referenceelements.hh>

#include <dune/fem/quadrature/caching/twistutility.hh>

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif // #if HAVE_UG

namespace Dune
{

  // UGGrid< 2 >
  // -----------

#if HAVE_UG
  template<>
  int TwistUtility< UGGrid< 2 > >::twistInSelf ( const GridType &grid, const LeafIntersection &it )
  {
    // for simplex twist is 0 
    // for cube twist is 1 for side 0 and 3 
    // for 1 and 2 is 0 
    if( it.inside()->type().isCube() )
    {
      const int face = it.indexInInside();
      return ((face == 1) || (face == 2) ? 0 : 1);
    }
    else
      return 0;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInSelf ( const LeafIntersection &it ) const
  {
    // for simplex twist is 0 
    // for cube twist is 1 for side 0 and 3 
    // for 1 and 2 is 0 
    if( it.inside()->type().isCube() )
    {
      const int face = it.indexInInside();
      return ((face == 1) || (face == 2) ? 0 : 1);
    }
    else
      return 0;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInSelf ( const GridType &grid, const LevelIntersection &it )
  {
    // for simplex twist is 0 
    // for cube twist is 1 for side 0 and 3 
    // for 1 and 2 is 0 
    if( it.inside()->type().isCube() )
    {
      const int face = it.indexInInside();
      return ((face == 1) || (face == 2) ? 0 : 1);
    }
    else
      return 0;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInSelf ( const LevelIntersection &it ) const
  {
    // for simplex twist is 0 
    // for cube twist is 1 for side 0 and 3 
    // for 1 and 2 is 0 
    if( it.inside()->type().isCube() )
    {
      const int face = it.indexInInside();
      return ((face == 1) || (face == 2) ? 0 : 1);
    }
    else
      return 0;
  }



  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const GridType &grid, const LeafIntersection &it )
  {
    assert( it.neighbor() );
    if( it.outside()->type().isCube() )
    {
      const int face = it.indexInOutside();
      return ((face == 1) || (face == 2) ? 1 : 0);
    }
    else
      return 1;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const LeafIntersection &it ) const
  {
    assert( it.neighbor() );
    if( it.outside()->type().isCube() )
    {
      const int face = it.indexInOutside();
      return ((face == 1) || (face == 2) ? 1 : 0);
    }
    else
      return 1;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const GridType &grid, const LevelIntersection &it )
  {
    assert( it.neighbor() );
    if( it.outside()->type().isCube() )
    {
      const int face = it.indexInOutside();
      return ((face == 1) || (face == 2) ? 1 : 0);
    }
    else
      return 1;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const LevelIntersection &it ) const
  {
    assert( it.neighbor() );
    if( it.outside()->type().isCube() )
    {
      const int face = it.indexInOutside();
      return ((face == 1) || (face == 2) ? 1 : 0);
    }
    else
      return 1;
  }



  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const GridType &grid, const LeafIntersection &it )
  {
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }

  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const LeafIntersection &it ) const
  { 
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }
    
  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const GridType &grid, const LevelIntersection &it )
  {
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }

  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const LevelIntersection &it ) const
  { 
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }
#endif // #if HAVE_UG



  // UGGrid< 3 >
  // -----------

#if HAVE_UG
  namespace UG3
  {

    struct CubeTwists 
    {
      static int twistInNeighbor ( const int face )
      {
        static const int twistInNeigh[6] = { 0, -2, -2,  0,  0, -1 };
        assert( face >= 0 && face < 6 );
        return twistInNeigh[ face ];
      }
      
      static int twistInSelf( const int face )
      {
        static const int twistInSelf[6]  = {-2,  0,  0, -2, -1,  0 };
        assert( face >= 0 && face < 6 );
        return twistInSelf[ face ];
      }
    };



    template< class IndexSet, class Entity >
    inline int calculateSimplexTwistInNeighbor
      ( const IndexSet &indexSet, const Entity &en, const int inEn,
        const Entity &nb, const int inNb )
    {
      typedef typename Entity::ctype ctype;
      static const int dim = Entity::dimension;

      const GenericReferenceElement< ctype, dim > &enRef
        = GenericReferenceElements< ctype, dim >::general( en.type() );
      const GenericReferenceElement< ctype, dim > &nbRef
        = GenericReferenceElements< ctype, dim >::general( nb.type() );
        
      // number of vertices of face 
      const int numVertices = enRef.size( inEn, 1, dim );
      assert( numVertices == nbRef.size( inNb, 1, dim ) );
      int enVx[ 4 ];
      int nbVx[ 4 ];
      
      int faceMap[ 4 ] = { 0, 1, 2, 3 };

      bool allRight = true;
      for( int i = 0; i < numVertices; ++i ) 
      {
        enVx[ i ] = indexSet.subIndex( en, enRef.subEntity( inEn, 1, i ,dim ), dim );
        nbVx[ i ] = indexSet.subIndex( nb, nbRef.subEntity( inNb, 1, i, dim ), dim );
        allRight &= (enVx[ i ] == nbVx[ i ]);
      }

      if( !allRight )
      {
        for( int i = 0; i < numVertices; ++i )
        {
          if( enVx[ i ] != nbVx[ i ] )
          {
            for( int k = 1; k < numVertices; ++k )
            {
              int newvx = (i+k) % numVertices;
              if( enVx[ i ] == nbVx[ newvx ] )
                faceMap[ i ] = newvx;
            }
          }
        }
      }

      const bool posOrientation = (faceMap[ 1 ] == (faceMap[ 0 ]+1) % numVertices);
      const int twist = posOrientation ? faceMap[ 0 ] : faceMap[ 1 ] - numVertices;

      const int mapTriangle[ 6 ] = {-2, -3, -1, 0, 2, 1};
      return (numVertices == 3 ? mapTriangle[ twist ] : twist);
    }

  }



  template<>
  int TwistUtility< UGGrid< 3 > >::twistInSelf ( const GridType &grid, const LeafIntersection &it )
  {
    if( it.inside()->type().isSimplex() )
    {
      // inside twist for simplices is zero 
      return 0;
    }
    else 
    {
      assert( it.inside()->type().isCube() );
      return UG3::CubeTwists::twistInSelf( it.indexInInside() );
    }
  }

  template<>
  int TwistUtility< UGGrid< 3 > >::twistInSelf ( const LeafIntersection &it ) const
  {
    DUNE_THROW( NotImplemented, "not implemented because grid is missing!" );
  }

  template<>
  int TwistUtility< UGGrid< 3 > >::twistInSelf ( const GridType &grid, const LevelIntersection &it )
  {
    if( it.inside()->type().isSimplex() )
    {
      // inside twist for simplices is zero 
      return 0;
    }
    else 
    {
      assert( it.inside()->type().isCube() );
      return UG3::CubeTwists::twistInSelf( it.indexInInside() );
    }
  }

  template<>
  int TwistUtility< UGGrid< 3 > >::twistInSelf ( const LevelIntersection &it ) const
  {
    DUNE_THROW( NotImplemented, "not implemented because grid is missing!" );
  }


  
  template<>
  int TwistUtility< UGGrid< 3 > >::twistInNeighbor ( const GridType &grid, const LeafIntersection &it )
  {
    assert( it.neighbor() );
    if( it.outside()->type().isSimplex() )
    {
      return UG3::calculateSimplexTwistInNeighbor
        ( grid.leafIndexSet(), *it.inside(), it.indexInInside(), *it. outside(), it.indexInOutside() );
    }
    else 
    {
      assert( it.outside()->type().isCube() );
      return UG3::CubeTwists::twistInNeighbor( it.indexInOutside() );
    }
  }

  template<>
  int TwistUtility< UGGrid< 3 > >::twistInNeighbor ( const LeafIntersection &it ) const
  {
    DUNE_THROW( NotImplemented, "not implemented because grid is missing!" );
  }

  template<>
  int TwistUtility< UGGrid< 3 > >::twistInNeighbor ( const GridType &grid, const LevelIntersection &it )
  {
    assert( it.neighbor() );
    if( it.outside()->type().isSimplex() )
    {
      return UG3::calculateSimplexTwistInNeighbor
        ( grid.leafIndexSet(), *it.inside(), it.indexInInside(), *it. outside(), it.indexInOutside() );
    }
    else 
    {
      assert( it.outside()->type().isCube() );
      return UG3::CubeTwists::twistInNeighbor( it.indexInOutside() );
    }
  }

  template<>
  int TwistUtility< UGGrid< 3 > >::twistInNeighbor ( const LevelIntersection &it ) const
  {
    DUNE_THROW( NotImplemented, "not implemented because grid is missing!" );
  }



  template<>
  bool TwistUtility< UGGrid< 3 > >::conforming ( const GridType &grid, const LeafIntersection &it )
  {
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }

  template<>
  bool TwistUtility< UGGrid< 3 > >::conforming ( const LeafIntersection &it ) const
  { 
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }
  
  template<>
  bool TwistUtility< UGGrid< 3 > >::conforming ( const GridType &grid, const LevelIntersection &it )
  {
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }

  template<>
  bool TwistUtility< UGGrid< 3 > >::conforming ( const LevelIntersection &it ) const
  { 
    return (!it.neighbor() || (it.inside()->level() == it.outside()->level())); 
  }
#endif // #if HAVE_UG

}
