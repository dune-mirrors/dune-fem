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
typedef Dune::Fem::FunctionSpace< GridType::ctype, double, dimw, dimw+2 > FuncSpace;
typedef Dune::Fem::DiscontinuousGalerkinSpace< FuncSpace, GridPartType, polOrd > DGSpaceType;
typedef Dune::Fem::DiscontinuousGalerkinSpace< FuncSpace, GridPartType, polOrd+1 > DGSpaceType2;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FuncSpace, GridPartType, polOrd > LagrangeSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< DGSpaceType > DiscreteFunctionType1;
typedef Dune::Fem::AdaptiveDiscreteFunction< DGSpaceType2 > DiscreteFunctionType2;
typedef Dune::Fem::AdaptiveDiscreteFunction< LagrangeSpaceType > DiscreteFunctionType3;

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

    std::vector<double> error;

    // do normal computations
    {
      Dune::GridPtr< GridType > gridptr( std::to_string( dimw ) + "dgrid.dgf" );
      GridType &grid = *gridptr;
      const int step = Dune::DGFGridInfo< GridType >::refineStepsForHalf();
      grid.globalRefine( step );

      GridPartType gridPart( grid );
      DGSpaceType space1( gridPart );
      DGSpaceType2 space2( gridPart );
      LagrangeSpaceType space3( gridPart );

      DiscreteFunctionType1 solution1( "sol1", space1 );
      DiscreteFunctionType2 solution2( "sol2", space2 );
      DiscreteFunctionType3 solution3( "sol3", space3 );

      // add solution to persistence manager to force backup
      Dune::Fem::persistenceManager << solution1;

      // TODO: fix multiple solution backup
      Dune::Fem::persistenceManager << solution2;
      Dune::Fem::persistenceManager << solution3;

      auto gridFunction = Dune::Fem::gridFunctionAdapter( f, gridPart, 5 );
      Dune::Fem::interpolate( gridFunction, solution1 );
      Dune::Fem::interpolate( gridFunction, solution2 );
      Dune::Fem::interpolate( gridFunction, solution3 );
      Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
      error.push_back( l2norm.distance( gridFunction, solution1 ) );
      error.push_back( l2norm.distance( gridFunction, solution2 ) );
      error.push_back( l2norm.distance( gridFunction, solution3 ) );

      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout << "Initial L2-error: ";
        for( size_t i=0; i<error.size() ; ++i )
          std::cout << error[i] << " ";
        std::cout << std::endl;

        std::cout << "Writing checkpoint ... "<<std::endl;
      }

      // write checkpoint
      double time = 0;
      Dune::Fem::CheckPointer< GridType >::writeSingleCheckPoint( grid, time, true );
    }

    // reset PersistenceManager to initial state (otherwise the restore will not work)
    Dune::Fem::persistenceManager.reset();

    // restore from checkpoint
    {
      std::string checkpointfile( DATA_PATH );
      checkpointfile += "/data/checkpoint";
      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout <<"Reading Grid from checkpoint ... "<< std::endl;
      }
      Dune::GridPtr< GridType > gridPtr( Dune::Fem::CheckPointer< GridType >::restoreGrid( checkpointfile ) );
      GridType &grid = *gridPtr;

      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout <<"Reading Data from checkpoint ... "<< std::endl;
      }

      GridPartType gridPart( grid );
      DGSpaceType space1( gridPart );
      DGSpaceType2 space2( gridPart );
      LagrangeSpaceType space3( gridPart );

      DiscreteFunctionType1 solution1( "sol1", space1 );
      DiscreteFunctionType2 solution2( "sol2", space2 );
      DiscreteFunctionType3 solution3( "sol3", space3 );

      // add solution to persistence manager to force backup
      Dune::Fem::persistenceManager << solution1;
      Dune::Fem::persistenceManager << solution2;
      Dune::Fem::persistenceManager << solution3;

      solution1.clear();
      solution2.clear();
      solution3.clear();

      Dune::Fem::CheckPointer< GridType >::restoreData( grid, checkpointfile );

      solution1.communicate();
      solution2.communicate();
      solution3.communicate();

      Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
      auto gridFunction = Dune::Fem::gridFunctionAdapter( f, gridPart, 5 );
      std::vector<double> errorNew;
      errorNew.push_back( l2norm.distance( gridFunction, solution1 ) );
      errorNew.push_back( l2norm.distance( gridFunction, solution2 ) );
      // TODO: Fix dofcompression and restore for LagrangeSpace
      // errorNew.push_back( l2norm.distance( gridFunction, solution3 ) );

      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout <<"Check L2-norm ... "<< std::endl;
      }
      for( size_t i=0; i<errorNew.size(); ++i )
      {
        const double diff = std::abs( error[i] - errorNew[i] );
        if( Dune::Fem::Parameter::verbose() )
        {
          std::cout << errorNew[i] << std::endl;
          std::cout << "Difference: " << std::scientific << diff << std::endl;
        }

        if( diff > 1e-12 )
          DUNE_THROW(Dune::IOError,"restored solution yields wrong L2-error");
      }

      if( Dune::Fem::MPIManager::rank() == 0 )
        std::cout << "Successful!"<<std::endl;
    }

    return 0;
  }
  catch( const Exception& e )
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
