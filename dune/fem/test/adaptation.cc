#ifdef ALBERTAGRID
#undef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#undef WORLDDIM
#define WORLDDIM ALBERTA_DIM
#endif

// only perform this test for the 3d version of ALUGrid
#if defined ALUGRID_CONFORM || defined ALUGRID_SIMPLEX || defined ALUGRID_CUBE || defined UGGRID
#define RUN_PROGRAM
#endif

#if defined ALUGRID_CONFORM || defined UGGRID
#define CONFORMING_SPACE
#undef USECOMBINEDSPACE
//#define USECOMBINEDSPACE
#endif

#include <config.h>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/space/padaptivespace.hh>

#ifdef HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "poisson-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
};

// Scheme
// ------

template < class GridPart, class DiscreteFunctionSpaceType >
struct Scheme
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridType GridType;

  typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

#ifdef HAVE_DUNE_ISTL
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#else
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#endif

  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >  RestrictionProlongationType;

  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;

  Scheme( GridPartType &gridPart, const int step  = 0 )
    : gridPart_( gridPart ),
      grid_( gridPart_.grid() ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      restrictProlong_( solution_ ),
      adaptationManager_( gridPart_.grid(), restrictProlong_ ),
      step_( step )
  {
    solution_.clear();
  }

  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  template <class GF>
  void initialize( const GF& gridFunction )
  {
    Dune::Fem::interpolate( gridFunction, solution_ );
  }

  //! mark elements for adaptation
  bool mark ( double time ) const
  {
    int marked = 0;
    int count = 0;
    int total = 0;

    // loop over all elements
    for( const auto& entity : discreteSpace_ )
    {
      // find center
      auto center = entity.geometry().center();
      DomainType x = DomainType(-time);
      x += center;

      const int maxLevel_ = 3 * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
      // refine if 0.3 < |x| < 1.0, otherwise (possibly) coarsen
      if( x.two_norm() > 0.3 && x.two_norm() < 0.4 && entity.level() <= maxLevel_ + step_ )
      {
        grid_.mark( 1, entity );
        marked = 1;
        count ++;
      }
      else
      {
        grid_.mark( -1, entity );
      }

      total++;
    }

    // get global max
    marked = grid_.comm().max( marked );

    // print info
    if( static_cast< bool >( marked ) )
      std::cout << "P" << Dune::Fem::MPIManager::rank() << ": " << "marked (" << count << " of " << total << ")" << std::endl;

    return static_cast< bool >(marked);
  }

  //! do the adaptation for a given marking
  void adapt()
  {
    // apply adaptation and load balancing
    adaptationManager_.adapt();
  }

protected:
  GridPartType &gridPart_; // grid part(view), e.g. here the leaf grid the discrete space is build with
  GridType &grid_;
  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown
  RestrictionProlongationType restrictProlong_ ; // local restriction/prolongation object
  AdaptationManagerType  adaptationManager_ ;    // adaptation manager handling adaptation
  const int step_;
};

template< class FunctionSpace >
struct Function : Dune::Fem::Function< FunctionSpace, Function< FunctionSpace > >
{
  void evaluate( const typename FunctionSpace::DomainType &x, typename FunctionSpace::RangeType &y ) const
  {
    y[ 0 ] = x[ 0 ] * (1-x[0]);
  }
};

// algorithm
// ---------

template <bool combinedSpace, class HGridType>
double algorithm ( HGridType &grid, const int step )
{
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;
  GridPartType gridPart(grid);

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, HGridType :: dimensionworld, 1 > SingleFunctionSpace;

  typedef Dune::Fem::DiscontinuousGalerkinSpace< SingleFunctionSpace, GridPartType, POLORDER > DGSpaceType;

  typedef typename
  std::conditional< combinedSpace,
      Dune::Fem::TupleDiscreteFunctionSpace< DGSpaceType, DGSpaceType >,
#ifdef CONFORMING_SPACE
      Dune::Fem::LagrangeDiscreteFunctionSpace< SingleFunctionSpace, GridPartType, POLORDER >
#else
      DGSpaceType
#endif
        >::type DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType  FunctionSpaceType;

  typedef Scheme< GridPartType, DiscreteFunctionSpaceType > SchemeType;
  SchemeType scheme( gridPart, step );

  typedef Function< FunctionSpaceType > FunctionType;
  FunctionType f;
  typedef Dune::Fem::GridFunctionAdapter< FunctionType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", f, gridPart, 5 );

  scheme.initialize( gridExactSolution );

  // output
  typedef std::tuple< const typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution);
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  for( double time = 0; time <= 0.15; time += 0.05 )
  {
    if( Dune::Fem::MPIManager::rank() == 0 )
      std::cout << "time: " << time << std::endl;

    // mark element for adaptation
    int max = 0;
    while( scheme.mark( time ) )
    {
      // adapt grid
      scheme.adapt();

      max++;
      if( max > 3 )
        break;
    }

    // data I/O
    dataOutput.write();
  }

  return 0.0;
}


int main ( int argc, char **argv )
try
{
#ifndef RUN_PROGRAM
  std::cerr << "No ALUGRID_CONFORM and GRIDDIM 3, so do nothing" << std::endl;
  return 0;
#endif

  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // this test only seems to work in single thread mode
  // not available at the moment: Dune::Fem::MPIManager::setMaxNumberThreads( 1 );
  Dune::Fem::MPIManager::setNumThreads( 1 );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );
  // make sure output parameters are added
  Dune::Fem::Parameter::append( "fem.prefix","output" );
  Dune::Fem::Parameter::append( "fem.io.savecount", "1" );
  Dune::Fem::Parameter::append( "fem.io.outputformat", "vtk-cell" );
  Dune::Fem::Parameter::append( "fem.io.partitioning", "rank" );
  Dune::Fem::Parameter::append( "fem.loadbalancing.step", "1" );
  Dune::Fem::Parameter::append( "fem.adaptation.method", "callback" );

  // type of hierarchical grid
  typedef Dune :: GridSelector :: GridType  HGridType ;

  // create grid from DGF file
  std::stringstream dgf;
  dgf << "DGF" << std::endl;
  dgf << "INTERVAL" << std::endl;
  dgf << "0 0 0" << std::endl;
  dgf << "1 1 1" << std::endl;
  const int nxcells = 2;
  dgf << nxcells << " " << nxcells << " " << nxcells << std::endl;
  dgf << "#" << std::endl;

  // the method rank and size from MPIManager are static
  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( dgf );
  HGridType& grid = *gridPtr ;

  // initial load balance
  grid.loadBalance();

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "level", 0 );
  const int repeats = Dune::Fem::Parameter::getValue< int >( "repeats", 2 );

  // number of global refinements to bisect grid width
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // normal space
  {
    // refine grid
    Dune :: Fem :: GlobalRefine::apply( grid, level * refineStepsForHalf );

    // calculate first step
    for( int step = 0; step < repeats; ++step )
      algorithm< false >( grid, step );
  }

  // combined space
  /*
  {
    // coarsen
    grid.globalRefine( -level*refineStepsForHalf );
    // refine grid
    Dune :: Fem :: GlobalRefine::apply( grid, level * refineStepsForHalf );

    // calculate first step
    for( int step = 0; step < repeats; ++step )
      algorithm< true >( grid, step );
  }
  */
  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
