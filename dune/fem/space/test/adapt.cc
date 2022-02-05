#include <config.h>

#include <iostream>

using namespace Dune;

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/padaptivespace/padaptivespace.hh>
#include <dune/fem/space/hpdg/orthogonal.hh>
#include <dune/fem/space/combinedspace.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/capabilities.hh>

#if HAVE_GRAPE && WANT_GRAPE && GRIDDIM > 1
#define USE_GRAPE 1
#else
#define USE_GRAPE 0
#endif

#if USE_GRAPE && GRIDDIM > 1
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/io/parameter.hh>

#include <dune/fem/test/testgrid.hh>

// polynom approximation order of quadratures,
// at least poolynom order of basis functions
static const int polOrd   = POLORDER;
static const int dimRange = 3 ;

using namespace Fem;

//***********************************************************************
/*! L2 Projection of a function f:
*/
//***********************************************************************

//! the index set we are using
typedef GridSelector::GridType MyGridType;
typedef DGAdaptiveLeafGridPart< MyGridType > GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
// typedef MatrixFunctionSpace < double , double, MyGridType::dimensionworld, 2,5 > FuncSpace;

#ifdef USE_COMBINED_SPACE
  typedef DiscontinuousGalerkinSpace< FunctionSpace < double , double, MyGridType::dimensionworld, 1 >,
                                      GridPartType, polOrd, CachingStorage> ContainedDiscreteFunctionSpaceType;
  #ifdef POINTBASED
    static const std::string usingSpaceName("Using CombinedSpace< dimRange, PointBased >");
    typedef CombinedSpace< ContainedDiscreteFunctionSpaceType, dimRange, PointBased > DiscreteFunctionSpaceType ;
  #else
    static const std::string usingSpaceName("Using CombinedSpace< dimRange, VariableBased >");
    typedef CombinedSpace< ContainedDiscreteFunctionSpaceType, dimRange, VariableBased > DiscreteFunctionSpaceType ;
  #endif
#else
  #ifdef USE_TUPLE_SPACE
    static const std::string usingSpaceName("Using TupleDiscreteFunctionSpace");
    typedef Dune::Fem::FunctionSpace< typename MyGridType::ctype, double, MyGridType::dimensionworld, dimRange+2 > FuncSpace1;
    typedef Dune::Fem::FunctionSpace< typename MyGridType::ctype, double, MyGridType::dimensionworld, dimRange > FuncSpace2;
    typedef Dune::Fem::DiscontinuousGalerkinSpace< FuncSpace1, GridPartType, polOrd+1, CachingStorage > DiscreteFunctionSpaceType1;
    typedef Dune::Fem::DiscontinuousGalerkinSpace< FuncSpace2, GridPartType, polOrd, CachingStorage > DiscreteFunctionSpaceType2;
    typedef Dune::Fem::TupleDiscreteFunctionSpace< DiscreteFunctionSpaceType1, DiscreteFunctionSpaceType2 > DiscreteFunctionSpaceType;
  #else
    static const std::string usingSpaceName("Using DiscontinuousGalerkinSpace< dimRange >");
typedef DiscontinuousGalerkinSpace< FunctionSpace < double , double, MyGridType::dimensionworld, dimRange >,
                                    GridPartType, polOrd, CachingStorage>  DiscreteFunctionSpaceType;
//typedef hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace < double , double, MyGridType::dimensionworld, dimRange >,
//                                    GridPartType, polOrd, true >  DiscreteFunctionSpaceType;
//typedef PAdaptiveDGSpace< FunctionSpace < double , double, MyGridType::dimensionworld, dimRange >,
//                                    GridPartType, polOrd, CachingStorage >  DiscreteFunctionSpaceType;
  #endif
#endif


////typedef CombinedSpace< ContainedDiscreteFunctionSpaceType, dimRange, PointBased > DiscreteFunctionSpaceType ;

typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

//! define the type of discrete function we are using , see

#if HAVE_PETSC
// PetscDiscreteFunction uses AdaptiveDiscreteFunction for dof prolongation and
// resttriction
typedef PetscDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#else
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#endif


//typedef Dune::Fem::ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteFunctionSpaceType, std::vector< double > > > DiscreteFunctionType;
//typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//typedef Dune::Fem::ReferenceBlockVector< FunctionSpaceType::RangeFieldType,
//                                         DiscreteFunctionSpaceType::localBlockSize > BlockVectorType;
//typedef Dune::Fem::BlockVectorDiscreteFunction< DiscreteFunctionSpaceType, BlockVectorType > DiscreteFunctionType;

typedef DofManager< MyGridType > DofManagerType;

typedef AdaptationManager< MyGridType, RestrictProlongDefault< DiscreteFunctionType > > AdaptationManagerType;

// ***********************************************************
// the exact solution to the problem for EOC calculation
struct ExactSolution
: public Fem::Function< FunctionSpaceType, ExactSolution >
{
  typedef FunctionSpaceType::RangeType RangeType;
  typedef FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef FunctionSpaceType::DomainType DomainType;

  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate ( const DomainType &x, RangeType &ret ) const
  {
    ret = 2.; // maximum of function is 2
    for( int i = 0; i < DomainType::dimension; ++i )
      ret *= sin( x[ i ]*(1.0 -x[ i ])*4.);
  }

  void evaluate ( const DomainType &x, RangeFieldType time, RangeType &ret ) const
  {
    evaluate( x, ret );
  }
};

// ********************************************************************
void adapt( MyGridType &grid, DiscreteFunctionType &solution, int step )
{
  const DiscreteFunctionType::DiscreteFunctionSpaceType &space = solution.space();
  RestrictProlongDefault<DiscreteFunctionType> rp(solution);
  rp.setFatherChildWeight(DGFGridInfo< MyGridType >::refineWeight());

  AdaptationManagerType adop(grid,rp);

  std::string message;

  int mark = 1;
  int count = std::abs(step);

  if(step < 0)
  {
    message += "Coarsening...";
    mark = -1;
  }
  else
    message += "Refining...";

  if( Parameter :: verbose() )
  {
    std::cout << "Grid leaf size:             " << grid.size( 0 ) << std::endl;
    std::cout << "AdaptiveLeafIndexSet.size:  " << space.indexSet().size( 0 ) << std::endl;
  }

  for( int i = 0; i < count; ++i )
  {
    for( const auto& entity : space)
      grid.mark( mark, entity );
    adop.adapt();
    if( Parameter::verbose() )
      std::cout << message << std::endl;
  }

  if( Parameter :: verbose() )
  {
    std::cout << "Grid leaf size:             " << grid.size( 0 ) << std::endl;
    std::cout << "AdaptiveLeafIndexSet.size:  " << space.indexSet().size( 0 ) << std::endl;
  }

}
// ********************************************************************
double algorithm ( MyGridType &grid, DiscreteFunctionType &solution, int step, int turn )
{
  const unsigned int order = solution.space().order();

  DiscreteFunctionType tmp ( solution );

  ExactSolution f;
  auto gridFunc = gridFunctionAdapter(f, solution.space().gridPart(), 2);
  {
    Dune::Fem::interpolate( gridFunc, solution );
    Dune :: Fem :: L2Norm< GridPartType > l2norm ( solution.space().gridPart(), 2*order+2 ) ;
    double new_error = l2norm.distance( f ,solution );
    if( Parameter::verbose() )
      std::cout << "before ref." << new_error << "\n\n";
  }

  adapt(grid,solution,step);

  // tmp solution should be zero after adapt
  double tmpErr = tmp.normSquaredDofs();
  if( tmpErr > 0 )
  {
    // return big error
    return 1e80;
  }

#if USE_GRAPE
  // if Grape was found, then display last solution
  if(turn > 0)
  {
    std::cerr << "GRAPE 1" << std::endl;
    GrapeDataDisplay < MyGridType > grape(grid);
    grape.dataDisplay( solution );
  }
#endif

  // calculation L2 error on refined grid
  // pol ord for calculation the error chould by higher than
  // pol for evaluation the basefunctions
  Dune :: Fem :: L2Norm< GridPartType > l2norm ( solution.space().gridPart(), 2*order+2 ) ;
  double error = l2norm.distance( gridFunc, solution );

#if USE_GRAPE
  // if Grape was found, then display last solution
  if(turn > 0)
  {
    std::cerr << "GRAPE 2" << std::endl;
    GrapeDataDisplay< MyGridType > grape(grid);
    grape.dataDisplay( solution );
  }
#endif

  //! perform l2-projection to refined grid
  Dune::Fem::interpolate( gridFunc, solution );
  double new_error = l2norm.distance( gridFunc, solution );
  if( Parameter::verbose() )
  {
    std::cout << "\nL2 Error : " << error << " on new grid " << new_error << "\n\n";
  }
#if USE_GRAPE
  // if Grape was found, then display last solution
  if(turn > 0)
  {
    std::cerr << "SIZE: " << solution.space().size()
        << " GRID: " << grid.size(0) << std::endl;
    std::cerr << "GRAPE 3" << std::endl;
    GrapeDataDisplay< MyGridType > grape(grid);
    grape.dataDisplay( solution );
  }
#endif
  return error;
}

//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main( int argc, char *argv[] )
try {
  MPIManager :: initialize( argc, argv );

  const char* paramName = "parameter";
  if( argc < 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << "<parameter>" << std :: endl;
  }
  else
    paramName = argv[1];

  std::string paramFile( paramName );
  std::cout << usingSpaceName << std::endl;

  // append parameter
  Parameter :: append( argc , argv );
  Parameter :: append( paramFile );

  int ml = 2 ; // default value = 2
  ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  std::vector<double> error(ml);

  MyGridType &grid = Dune::Fem::TestGrid::grid();
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();

  GridPartType part ( grid );
  DiscreteFunctionSpaceType space( part );

  // threshold for EOC difference to predicted value
  const double eocThreshold = Parameter :: getValue("adapt.eocthreshold", double(0.2) );

  const bool isLocallyAdaptive = Dune::Fem::Capabilities::isLocallyAdaptive< MyGridType > :: v ;

  DiscreteFunctionType solution ( "sol", space );
  solution.clear();
  if( Parameter::verbose() )
    std::cout << "------------    Refining:" << std::endl;
  for(int i=0; i<ml; i+=1)
  {
    error[i] = algorithm ( grid , solution, step, (i==ml-1));
    if (i>0)
    {
      if ( isLocallyAdaptive )
      {
        double eoc = log( error[i-1]/error[i]) / M_LN2;
        if( Parameter::verbose() )
          std::cout << "EOC = " << eoc << std::endl;
        if( std::abs( eoc - (space.order()+eocThreshold) ) < 0 )
        {
          DUNE_THROW(InvalidStateException,"EOC check of refinement failed");
        }
      }
      else
        std::cout << "no EOC for non-adaptive grid" << std::endl;
    }
  }
  if( Parameter::verbose() )
    std::cout << "------------   Coarsening:" << std::endl;
  for(int i=ml-1; i>=0; i-=1)
  {
    error[i] = algorithm ( grid , solution,-step, 1);
    if (i<ml-1)
    {
      if( isLocallyAdaptive )
      {
        double eoc = log( error[i+1]/error[i]) / M_LN2;
        if( Parameter::verbose() )
          std::cout << "EOC = " << eoc << std::endl;
        if( std::abs( eoc + (space.order()+eocThreshold) ) < 0 )
        {
          DUNE_THROW(InvalidStateException,"EOC check of coarsening failed");
        }
      }
      else
        std::cout << "no EOC for non-adaptive grid" << std::endl;
    }
  }
  return 0;
}
catch( const Dune :: Exception &exception )
{
  std :: cerr << exception << std :: endl;
  return 1;
}
