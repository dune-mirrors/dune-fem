#include "basefunctiontest.hh"

#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/quadrature/cachequad.hh>

namespace Dune {

  void LagrangeBase_Test::run() 
  {
    testBaseFunctions();
  }

  void LagrangeBase_Test::testBaseFunctions() 
  {
    GridPtr<GridType> gridPtr(gridFile_); 
    GridType& grid = *gridPtr;

    typedef LeafGridPart<GridType> GridPartType; 
    LeafGridPart<GridType> gridPart(grid);

    // check polynomial order 1 
    typedef FunctionSpace<double,double,dimworld,1> FunctionSpaceType;
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType,
            GridPartType, 1 > OneSpaceType;

    {
      OneSpaceType space(gridPart);
      checkLagrangeBase( space );
    }
    
    /*
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType,
            GridPartType, 2 > TwoSpaceType;

    {
      TwoSpaceType space(gridPart);
      checkLagrangeBase( space );
    }
    */
  }
  
  template <class SpaceType> 
  void LagrangeBase_Test::checkLagrangeBase(const SpaceType& space) 
  {
    typedef typename SpaceType :: IteratorType IteratorType; 
    typedef typename SpaceType :: BaseFunctionSetType
      BaseFunctionSetType;
    typedef typename SpaceType :: DomainType DomainType; 
    typedef typename SpaceType :: RangeType RangeType; 
    typedef typename SpaceType :: GridType :: template
      Codim<0>::Geometry Geometry; 
    
    IteratorType end = space.end();
    for(IteratorType it = space.begin(); it != end; ++it)
    {
      const BaseFunctionSetType& baseSet = space.baseFunctionSet( *it );
      const int numBaseFct = baseSet.numBaseFunctions();

      const Geometry& geo = it->geometry();

      for(int i=0; i<numBaseFct; ++i) 
      {
        RangeType phi(0.0);     
        assert( geo.corners() == numBaseFct );
        DomainType point = geo.local( geo[i] );
        
        // eval on lagrange point 
        baseSet.eval( i , point , phi ); 
        //baseSet.eval( i , quad, i, phi ); 

        assert( std::abs(phi[0] - 1.0) < 1e-10 );
        
        for(int j=0; j<numBaseFct; ++j) 
        {
          if(i == j) continue;

          // eval on lagrange point 
          baseSet.eval( j , point , phi ); 
          //baseSet.eval( i , quad, j, phi ); 

          assert( std::abs(phi[0]) < 1e-10 );
        }
      }
    }
  }

} // end namespace Dune
