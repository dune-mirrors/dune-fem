#ifdef ALBERTAGRID
// set dimensions to ALBERTA dimensions to avoid conflicts
#undef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#undef WORLDDIM
#define WORLDDIM ALBERTA_DIM
#endif

#ifdef YASPGRID
#define SKIP_TEST_BECAUSE_USING_YASPGRID
#endif

#include <config.h>

#include <iostream>

using namespace Dune;

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/l2norm.hh>

#if HAVE_GRAPE && GRIDDIM > 1
#define USE_GRAPE 1
#else
#define USE_GRAPE 0
#endif

#if USE_GRAPE && GRIDDIM > 1
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/io/parameter.hh>

#include <dune/fem/test/testgrid.hh>

using namespace Fem;

// polynom approximation order of quadratures,
// at least poolynom order of basis functions
const int polOrd = POLORDER;

#ifndef GRIDDIM
#define GRIDDIM dimworld
#endif

//***********************************************************************
/*! L2 Projection of a function f:
*/
//***********************************************************************

//! the index set we are using
typedef GridSelector::GridType MyGridType;
typedef DGAdaptiveLeafGridPart< MyGridType , All_Partition > GridPartType;
//typedef AdaptiveLeafGridPart< MyGridType > GridPartType;
//typedef HierarchicGridPart< MyGridType > GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
// typedef MatrixFunctionSpace < double , double, GRIDDIM , 2,5 > FuncSpace;

//! define the function space our unkown belong to
//! see dune/fem/lagrangebase.hh
typedef FunctionSpace < double , double, GRIDDIM , 5 > FuncSpace;
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType,
  polOrd,CachingStorage> DiscreteFunctionSpaceType;

typedef FunctionSpace < double , double, GRIDDIM , 3 > FuncSpace2;
typedef DiscontinuousGalerkinSpace<FuncSpace2, GridPartType,
  polOrd,CachingStorage> DiscreteFunctionSpaceType2;

typedef LegendreDiscontinuousGalerkinSpace<FuncSpace2, GridPartType,
  polOrd,CachingStorage> DiscreteFunctionSpaceType3;

//! define the type of discrete function we are using , see
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef DofManager< MyGridType > DofManagerType;

typedef AdaptationManager< MyGridType, RestrictProlongDefault< DiscreteFunctionType > > AdaptationManagerType;

// ***********************************************************
// the exact solution to the problem for EOC calculation
struct ExactSolution
: public Fem::Function< FuncSpace, ExactSolution >
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;

  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate ( const DomainType &x, RangeType &ret ) const
  {
    ret = 2.; // maximum of function is 2
    for( int i = 0; i < DomainType::dimension; ++i )
      ret *= sin( x[ i ]*(1.0 -x[ i ])*4.);
  }
};

// ********************************************************************
class DGL2ProjectionAllPartitionNoComm
{

 public:
  template <class FunctionType, class DiscreteFunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc,
                       int polOrd = -1 )
  {
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    const DiscreteFunctionSpaceType& space =  discFunc.space();

    if( polOrd < 0 )  polOrd = 2 * space.order() + 2 ;

    discFunc.clear();

    TemporaryLocalFunction< DiscreteFunctionSpaceType > dfLocal( space );

    typename DiscreteFunctionSpaceType::RangeType ret (0.0);
    typename DiscreteFunctionSpaceType::RangeType phi (0.0);

    for( const auto& entity : elements( space.gridPart(), Partitions::all ) )
    {
      const auto &itGeom = entity.geometry();
      dfLocal.bind( entity );
      dfLocal.clear();

      CachingQuadrature<GridPartType,0> quad(entity, polOrd);
      for( size_t qP = 0; qP < quad.nop(); ++qP )
      {
        f.evaluate(itGeom.global(quad.point(qP)), ret);
        ret *= quad.weight(qP) ;
        dfLocal.axpy( quad[qP], ret );
      }

      discFunc.setLocalDofs(entity, dfLocal );
      dfLocal.unbind();
    }
  }
};

void resetNonInterior( DiscreteFunctionType &solution )
{
  typedef DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  const DiscreteFunctionSpaceType& space = solution.space();
  TemporaryLocalFunction< DiscreteFunctionSpaceType > dfLocal( space );

  int count = 0;
  for( const auto& entity : elements( space.gridPart(), Partitions::all ) )
  {
    if( entity.partitionType() != InteriorEntity )
    {
      ++count ;
      dfLocal.bind( entity );
      dfLocal.clear();
      solution.setLocalDofs( entity, dfLocal );
      dfLocal.unbind();
    }
  }

  std::cout << "P[" << space.gridPart().comm().rank() << "]  reset " << count << " entities "  << std::endl;

}

// ********************************************************************
double algorithm ( MyGridType &grid, DiscreteFunctionType &solution, int step, int turn )
{
  ExactSolution f;

  // create L2 Norm, no communication
  L2Norm< typename DiscreteFunctionType::GridPartType > l2norm( solution.space().gridPart(), 2*solution.space().order(), false );
  solution.clear();

  DGL2ProjectionAllPartitionNoComm :: project( f, solution );
  // compute l2 error on all elements
  double new_error = l2norm.distance( f ,solution, Partitions::all );
  std::cout << "P[" << grid.comm().rank() << "]  start comm: " << new_error << std::endl;

  // reset all non-interior entities,
  // these should be restored during communication
  resetNonInterior( solution );

  // do communication
  solution.communicate();

  // calculate l2 error again on all elements
  double error = l2norm.distance( f ,solution, Partitions::all );
  std::cout << "P[" << grid.comm().rank() << "]  done comm: " << error << std::endl;

  if( std::abs( new_error - error ) > 1e-10 )
    DUNE_THROW(InvalidStateException,"Communication not working correctly");

  ///////////////////////////////////////////////////
  //  test non-blocking communication
  ///////////////////////////////////////////////////

  typedef DiscreteFunctionType :: DiscreteFunctionSpaceType :: CommunicationManagerType
    :: NonBlockingCommunicationType  NonBlockingCommunicationType;

  NonBlockingCommunicationType nonBlocking =
    solution.space().communicator().nonBlockingCommunication();

  // send data
  nonBlocking.send( solution );

  // do some work,
  // reset all non-interior entities,
  // these should be restored during communication
  resetNonInterior( solution );

  // receive data
  nonBlocking.receive( solution );

  // calculate l2 error again
  double nonBlock = l2norm.distance(f ,solution, Partitions::all );
  std::cout << "P[" << grid.comm().rank() << "]  non-blocking: " << nonBlock << std::endl;

  if( std::abs( new_error - nonBlock ) > 1e-10 )
    DUNE_THROW(InvalidStateException,"Communication not working correctly");

  return error;
}

std::shared_ptr< GridPartType > gpPtr;
std::shared_ptr< DiscreteFunctionSpaceType > spcPtr;

std::shared_ptr< DiscreteFunctionType > createSolution(MyGridType& grid)
{
  {
    std::shared_ptr< GridPartType > gp1( new GridPartType(grid) );
    std::shared_ptr< DiscreteFunctionSpaceType > spc1( new DiscreteFunctionSpaceType(*gp1) );
    std::shared_ptr< DiscreteFunctionType > sol1( new DiscreteFunctionType("sol1", *spc1) );

    sol1->communicate();

    std::shared_ptr< GridPartType > gp2( new GridPartType(grid) );
    std::shared_ptr< DiscreteFunctionSpaceType > spc2( new DiscreteFunctionSpaceType(*gp2) );
    std::shared_ptr< DiscreteFunctionType > sol2( new DiscreteFunctionType("sol2", *spc2) );

    sol2->communicate();

    std::shared_ptr< DiscreteFunctionSpaceType2 > spc3( new DiscreteFunctionSpaceType2(*gp2) );
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType2 > DiscreteFunctionType2;
    std::shared_ptr< DiscreteFunctionType2 > sol3( new DiscreteFunctionType2("sol3", *spc3) );

    sol3->communicate();
    sol2->communicate();
  }

  gpPtr.reset( new GridPartType( grid ) );
  spcPtr.reset( new DiscreteFunctionSpaceType( *gpPtr ) );
  std::shared_ptr< DiscreteFunctionType > sol( new DiscreteFunctionType("sol", *spcPtr));
  return sol;
}

//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main( int argc, char *argv[] )
try {
#ifdef SKIP_TEST_BECAUSE_USING_YASPGRID
  std::cerr << "No comm check for YASPGRID because overlap not working correctly!" << std::endl;
  return 0;
#endif

  MPIManager :: initialize( argc, argv );

  const char* paramName = "parameter";
  if( argc < 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << "<parameter>" << std :: endl;
  }
  else
    paramName = argv[1];

  std::string paramFile( paramName );

  // append parameter
  Parameter :: append( argc , argv );
  Parameter :: append( paramFile );

  int ml = 2 ; // default value = 2
  //ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  std::vector<double> error(ml);

  MyGridType &grid = Dune::Fem::TestGrid::grid();
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();

  std::shared_ptr< DiscreteFunctionType > solution = createSolution( grid );

  for(int i=0; i<ml; ++i )
  {
    if( grid.comm().rank() == 0)
      std::cout << std::endl << "**** Start communication cycle " << i << "  ****" << std::endl;

    GlobalRefine :: apply( grid, step );
    error[ i ] = algorithm ( grid , *solution, step, (i==ml-1));
  }

  solution.reset();
  spcPtr.reset();
  gpPtr.reset();

  return 0;
}
catch( const Dune :: Exception &exception )
{
  std :: cerr << exception << std :: endl;
  return 1;
}

