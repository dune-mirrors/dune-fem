#include <config.h>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

// #include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
static const int dimw = Dune::GridSelector::dimworld;

#include <dune/common/exceptions.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>

using namespace Dune;
using namespace Fem;


// polynom approximation order of quadratures,
// at least polynom order of basis functions
const int polOrd = POLORDER;

//***********************************************************************
/*! L2 Projection of a function f:

   This is an example how to solve the equation on
   \f[\Omega = (0,1)^2 \f]

   \f[ \int_{\Omega} u \phi = \int_{\Omega} f \phi  \ \ \ in \Omega \f]
   \f[ f(x,y) = x ( 1 - x) y ( 1 - y ) \f]

   Here u is the L_2 projection of f.

   The Projection should converge to the given function f.
   with the finite element method using lagrangian elements of polynom order +1.
 */
//***********************************************************************
typedef Dune::GridSelector::GridType GridType;

//typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
typedef Dune::Fem::FunctionSpace< GridType::ctype, double, dimw, 4 > FuncSpace;
typedef Dune::Fem::DiscontinuousGalerkinSpace< FuncSpace, GridPartType, polOrd > DiscreteFunctionSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

// the exact solution to the problem for EOC calculation
struct ExactSolution
{
  typedef FuncSpace FunctionSpaceType;
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;

  ExactSolution () {}

  void evaluate ( const DomainType &x, RangeType &ret )  const
  {
    ret = 2.;
    for( int i = 0; i < DomainType::dimension; i++ )
      ret[ 0 ] *= x[ i ]*(1.0 -x[ i ]);
    for( int i = 0; i < DomainType::dimension; i++ )
      ret[ 3 ] *= sin( 2.*M_PI*x[ i ]*(1.0 -x[ i ]));
    ret[ 1 ] = -x[ 1 ];
    ret[ 2 ] = x[ 0 ];
  }
};



//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  Dune::Fem::Parameter::append( argc, argv );
  try
  {
    Dune::Fem::Parameter::append( "parameter" );
    ExactSolution f;

    double error = 0.0;

    // typedef std::tuple< DiscreteFunctionType * > DataTuple;

    // do normal computations
    {
      Dune::GridPtr< GridType > gridptr( std::to_string( dimw ) + "dgrid.dgf" );
      GridType &grid = *gridptr;
      const int step = Dune::DGFGridInfo< GridType >::refineStepsForHalf();
      grid.globalRefine( step );

      GridPartType gridPart( grid );
      DiscreteFunctionSpaceType space( gridPart );
      DiscreteFunctionType solution( "sol", space );

      // add solution to persistence manager to force backup
      Dune::Fem::persistenceManager << solution;

      auto gridFunction = Dune::Fem::gridFunctionAdapter( f, gridPart, 14 );
      Dune::Fem::interpolate( gridFunction, solution );
      Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
      error = l2norm.distance( gridFunction, solution );

      std::cout << "Initial L2-error: " << error << std::endl;

      std::cout << "Writing checkpoint ... "<<std::endl;
      // write checkpoint
      double time = 0;
      Dune::Fem::CheckPointer< GridType >::writeSingleCheckPoint( grid, time, true );
    }

    // reset PersistenceManager to initial state (otherwise the restore will not work)
    Dune::Fem::persistenceManager.reset();

    // restore from checkpoint
    {
      std::cout <<"Reading Grid from checkpoint ... "<< std::endl;
      std::string checkpointfile = "data/checkpoint";
      Dune::GridPtr< GridType > gridPtr( Dune::Fem::CheckPointer< GridType >::restoreGrid( checkpointfile ) );
      GridType &grid = *gridPtr;

      std::cout <<"Reading Data from checkpoint ... "<< std::endl;
      GridPartType gridPart( grid );
      DiscreteFunctionSpaceType space( gridPart );
      DiscreteFunctionType solution ( "sol", space );
      solution.clear();

      // register solution to persistence manager
      Dune::Fem::persistenceManager << solution;

      Dune::Fem::CheckPointer< GridType >::restoreData( grid, checkpointfile );

      std::cout <<"Check L2-norm ... "<< std::endl;
      Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
      auto gridFunction = Dune::Fem::gridFunctionAdapter( f, gridPart, 14 );
      const double errorNew =  l2norm.distance( gridFunction, solution );
      const double diff = std::abs( error - errorNew );
      std::cout << errorNew << std::endl;
      std::cout << "Difference: " << std::scientific << diff << std::endl;
      if( diff > 1e-12 )
        DUNE_THROW(Dune::IOError,"restored solution yields wrong L2-error");

      std::cout << "Successful!"<<std::endl;
    }

    return 0;
  }
  catch( Exception e )
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
