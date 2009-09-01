#include <config.h>

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
    return (it.inside()->type().isSimplex()) ? 0 : 
      (it.numberInSelf() == 1 || it.numberInSelf() == 2) ? 0 : 1;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInSelf ( const LeafIntersection &it ) const
  {
    // for simplex twist is 0 
    // for cube twist is 1 for side 0 and 3 
    // for 1 and 2 is 0 
    return (it.inside()->type().isSimplex()) ? 0 : 
      (it.numberInSelf() == 1 || it.numberInSelf() == 2) ? 0 : 1;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInSelf ( const GridType &grid, const LevelIntersection &it )
  {
    // for simplex twist is 0 
    // for cube twist is 1 for side 0 and 3 
    // for 1 and 2 is 0 
    return (it.inside()->type().isSimplex()) ? 0 : 
      (it.numberInSelf() == 1 || it.numberInSelf() == 2) ? 0 : 1;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInSelf ( const LevelIntersection &it ) const
  {
    // for simplex twist is 0 
    // for cube twist is 1 for side 0 and 3 
    // for 1 and 2 is 0 
    return (it.inside()->type().isSimplex()) ? 0 : 
      (it.numberInSelf() == 1 || it.numberInSelf() == 2) ? 0 : 1;
  }



  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const GridType &grid, const LeafIntersection &it )
  {
    assert( it.neighbor() );
    return (it.outside()->type().isSimplex()) ? 1 : 
      (it.numberInNeighbor() == 1 || it.numberInNeighbor() == 2) ? 1 : 0;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const LeafIntersection &it ) const
  {
    assert( it.neighbor() );
    return (it.outside()->type().isSimplex()) ? 1 : 
      (it.numberInNeighbor() == 1 || it.numberInNeighbor() == 2) ? 1 : 0;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const GridType &grid, const LevelIntersection &it )
  {
    assert( it.neighbor() );
    return (it.outside()->type().isSimplex()) ? 1 : 
      (it.numberInNeighbor() == 1 || it.numberInNeighbor() == 2) ? 1 : 0;
  }

  template<>
  int TwistUtility< UGGrid< 2 > >::twistInNeighbor ( const LevelIntersection &it ) const
  {
    assert( it.neighbor() );
    return (it.outside()->type().isSimplex()) ? 1 : 
      (it.numberInNeighbor() == 1 || it.numberInNeighbor() == 2) ? 1 : 0;
  }



  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const GridType &grid, const LeafIntersection &it )
  {
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true; 
  }

  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const LeafIntersection &it ) const
  { 
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true; 
  }
    
  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const GridType &grid, const LevelIntersection &it )
  {
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true; 
  }

  template<>
  bool TwistUtility< UGGrid< 2 > >::conforming ( const LevelIntersection &it ) const
  { 
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true; 
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
        return twistInNeigh[face];
      }
      
      static int twistInSelf( const int face )
      {
        static const int twistInSelf[6]  = {-2,  0,  0, -2, -1,  0 };
        assert( face >= 0 && face < 6 );
        return twistInSelf[face];
      }
    };



    template< class IndexSetType, class EntityType >
    inline int calculateSimplexTwistInNeighbor
      ( const IndexSetType &set, const EntityType &en, const int inSelf,
        const EntityType &nb, const int inNeigh )
    {
      typedef typename EntityType::ctype ctype;
      static const int dim = EntityType::dimension;

      const ReferenceElement< ctype, dim > &enRef
        = ReferenceElements< ctype, dim >::general( en.type() );
      const ReferenceElement< ctype, dim > &nbRef
        = ReferenceElements< ctype, dim >::general( nb.type() );
        
      // number of vertices of face 
      const int numVert = enRef.size (inSelf,1, dim);
      int enVx [4];
      int nbVx [4];
      
      int faceMap[4] = { 0, 1, 2, 3};

      bool allRight = true;
      for(int i=0; i<numVert; ++i ) 
      {
        enVx[i] = set.template subIndex<dim>(en, enRef.subEntity(inSelf,1,i,dim)); 
        nbVx[i] = set.template subIndex<dim>(nb, nbRef.subEntity(inNeigh,1,i,dim)); 
        if( enVx[i] != nbVx[i] ) allRight = false;
      }

      if( !allRight )
      {
        for(int i=0; i<numVert; ++i)
        {
          if(enVx[i] != nbVx[i])
          {
            for(int k=1; k<numVert; ++k)
            {
              int newvx = (i+k) % numVert;
              if( enVx[i] == nbVx[newvx] ) faceMap[i] = newvx;
            }
          }
        }
      }

      // return twist 

      if (faceMap[1] == (faceMap[0]+1) % numVert) 
      {
        return faceMap[0];
      }
      else 
      {
        int twst = faceMap[1] - numVert;
        if( numVert == 3 ) 
        {
          // same bug as in Alberta (check reference elements)
          if( twst == -3 ) return -2;
          else if ( twst == -2 ) return -3; 
          else return twst; 
        }
        else 
          return twst; 
      }
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
      return UG3::CubeTwists::twistInSelf( it.numberInSelf() );
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
      return UG3::CubeTwists::twistInSelf( it.numberInSelf() );
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
        ( grid.leafIndexSet(), *it.inside(), it.numberInSelf(), *it. outside(), it.numberInNeighbor() );
    }
    else 
    {
      assert( it.outside()->type().isCube() );
      return UG3::CubeTwists::twistInNeighbor( it.numberInNeighbor() );
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
        ( grid.leafIndexSet(), *it.inside(), it.numberInSelf(), *it. outside(), it.numberInNeighbor() );
    }
    else 
    {
      assert( it.outside()->type().isCube() );
      return UG3::CubeTwists::twistInNeighbor( it.numberInNeighbor() );
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
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true;
  }

  template<>
  bool TwistUtility< UGGrid< 3 > >::conforming ( const LeafIntersection &it ) const
  { 
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true;
  }
  
  template<>
  bool TwistUtility< UGGrid< 3 > >::conforming ( const GridType &grid, const LevelIntersection &it )
  {
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true;
  }

  template<>
  bool TwistUtility< UGGrid< 3 > >::conforming ( const LevelIntersection &it ) const
  { 
    return (it.neighbor()) ? 
      (it.inside()->level() == it.outside()->level()) : true;
  }
#endif // #if HAVE_UG

}
