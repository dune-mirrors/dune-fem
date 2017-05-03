#include "basisfunctiontest.hh"

// #include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/space/lagrange.hh>

namespace Dune
{

  namespace Fem
  {

    void LagrangeBasis_Test::run()
    {
      testBasisFunctions();
    }

    void LagrangeBasis_Test::testBasisFunctions()
    {
      typedef GridSelector::GridType GridType;
      static const int dimworld = GridSelector::dimworld;

      GridPtr< GridType > gridPtr( gridFile_ );
      GridType& grid = *gridPtr;

      typedef LeafGridPart< GridType > GridPartType;
      GridPartType gridPart( grid );

      typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;

      // check polynomial order 1
      typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >
        OneSpaceType;
      {
        std :: cout << "Linear Basis Functions" << std :: endl;
        OneSpaceType space( gridPart);
        checkLagrangeBasis( space );
      }

      #ifdef TEST_SECOND_ORDER
      // check polynomial order 2
      typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2 >
        TwoSpaceType;
      {
        std :: cout << "Quadratic Basis Functions" << std :: endl;
        TwoSpaceType space( gridPart );
        checkLagrangeBasis( space );
      }
      #endif

      #ifdef TEST_THIRD_ORDER
      typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 3 >
        ThreeSpaceType;
      {
        std :: cout << "Cubic Basis Functions" << std :: endl;
        ThreeSpaceType space( gridPart );
        checkLagrangeBasis( space );
      }
      #endif

    }

    template< class SpaceType >
    void LagrangeBasis_Test :: checkLagrangeBasis( const SpaceType &space )
    {
      typedef typename SpaceType :: BasisFunctionSetType  BasisFunctionSetType;
      typedef typename SpaceType :: DomainType            DomainType;
      typedef typename SpaceType :: RangeType             RangeType;
      typedef typename SpaceType :: GridPartType          GridPartType;

      typedef LagrangePointSet< GridPartType, SpaceType :: polynomialOrder > LagrangePointSetType;

      int errors = 0;
      for(const auto& entity : space)
      {
        const BasisFunctionSetType& basisSet = space.basisFunctionSet( entity );
        const int numBasisFct = basisSet.size();

        LagrangePointSetType pointSet( entity.type(), space.order( entity ) );
        const int numPoints = pointSet.size();

        std::vector< RangeType > phi( numPoints, RangeType( 0 ) );

        for( int i = 0; i < numPoints; ++i )
        {

          const DomainType& x = pointSet.point( i );

          // evaluate all basisFunctions on lagrange point i
          basisSet.evaluateAll( x, phi );

          if( std :: abs( phi[ i ][ 0 ] - 1.0 ) >= 1e-10 )
          {
            std :: cout << "Basis function " << i << " failed at " << x
                        << " (" << phi[ i ][ 0 ] << " != 1)!" << std :: endl;
            ++errors;
          }

          for( int j = 0; j < numBasisFct; ++j )
          {
            if( i == j )
              continue;

            // evaluate on lagrange point
            if( std :: abs( phi[ j ][ 0 ] ) >= 1e-10 )
            {
              std :: cout << "Basis function " << j << " failed at " << x
                          << " (" << phi[ j ][ 0 ] << " != 0)!" << std :: endl;
              ++errors;
            }
          }
        }
      }

      assert( errors == 0 );
    }

  } // namespace Fem

} // namespace Dune
