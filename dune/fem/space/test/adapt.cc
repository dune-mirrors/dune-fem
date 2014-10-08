#include <config.h>

#include <iostream>

using namespace Dune;

#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/function/adaptivefunction.hh>
//#include <dune/fem/function/vectorfunction.hh>
//#include <dune/fem/function/attachedfunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/operator/projection/dgl2projection.hh>
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

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
const int polOrd = POLORDER;

#ifndef GRIDDIM 
#define GRIDDIM dimworld 
#endif

using namespace Fem;

//***********************************************************************
/*! L2 Projection of a function f: 
*/
//***********************************************************************

//! the index set we are using
typedef GridSelector::GridType MyGridType;
typedef DGAdaptiveLeafGridPart< MyGridType > GridPartType;
//typedef AdaptiveLeafGridPart< MyGridType > GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
// typedef MatrixFunctionSpace < double , double, GRIDDIM , 2,5 > FuncSpace;

//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh
typedef FunctionSpace < double , double, GRIDDIM , 5 > FuncSpace;
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
                  polOrd,CachingStorage> DiscreteFunctionSpaceType;

//typedef LegendreDiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
//                  polOrd,CachingStorage> DiscreteFunctionSpaceType;

//typedef LagrangeDiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
//                  polOrd,CachingStorage> DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

//typedef ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteFunctionSpaceType, DynamicVector< double > > > DiscreteFunctionType;
//typedef AttachedDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

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

  void evaluate ( const DomainType &x, RangeFieldType time, RangeType &ret ) const
  {
    evaluate( x, ret );
  }
};
 
// ********************************************************************
void adapt( MyGridType &grid, DiscreteFunctionType &solution, int step )
{
  typedef DiscreteFunctionType::DiscreteFunctionSpaceType::IteratorType Iterator;
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
    const Iterator end = space.end();
    for( Iterator it = space.begin(); it != end ; ++it )
      grid.mark( mark, *it );
    adop.adapt();
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

  {
    ExactSolution f;
    DGL2ProjectionImpl :: project( f, solution );
    Dune :: Fem :: L2Norm< GridPartType > l2norm ( solution.space().gridPart(), 2*order+2 ) ;
    double new_error = l2norm.distance( f ,solution ); 
    std::cout << "before ref." << new_error << "\n\n"; 
  }

  adapt(grid,solution,step);

#if USE_GRAPE
  // if Grape was found, then display last solution 
  if(turn > 0) 
  {
    std::cerr << "GRAPE 1" << std::endl;
    GrapeDataDisplay < MyGridType > grape(grid); 
    grape.dataDisplay( solution );
  }
#endif

  ExactSolution f; 
  // calculation L2 error on refined grid
  // pol ord for calculation the error chould by higher than 
  // pol for evaluation the basefunctions 
  Dune :: Fem :: L2Norm< GridPartType > l2norm ( solution.space().gridPart(), 2*order+2 ) ;
  double error = l2norm.distance( f, solution );

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
  DGL2ProjectionImpl :: project ( f, solution );
  double new_error = l2norm.distance( f, solution );
  std::cout << "\nL2 Error : " << error << " on new grid " << new_error << "\n\n";
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

  // append parameter 
  Parameter :: append( argc , argv );
  Parameter :: append( paramFile );

  int ml = 2 ; // default value = 2 
  ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  std::vector<double> error(ml);

  std::ostringstream gridFilenameStream;
  gridFilenameStream << GRIDDIM << "dgrid.dgf";
  GridPtr< MyGridType > grid( gridFilenameStream.str() );

  const int step = DGFGridInfo< MyGridType >::refineStepsForHalf();

  GridPartType part ( *grid );
  DiscreteFunctionSpaceType space( part );

  // threshold for EOC difference to predicted value 
  const double eocThreshold = Parameter :: getValue("adapt.eocthreshold", double(0.2) );

  const bool isLocallyAdaptive = Dune::Fem::Capabilities::isLocallyAdaptive< MyGridType > :: v ;

  DiscreteFunctionType solution ( "sol", space );
  solution.clear();
  std::cout << "------------    Refining:" << std::endl;
  for(int i=0; i<ml; i+=1)
  {
    error[i] = algorithm ( *grid , solution, step, (i==ml-1));
    if (i>0) 
    {
      if ( isLocallyAdaptive ) 
      {
        double eoc = log( error[i-1]/error[i]) / M_LN2; 
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
  std::cout << "------------   Coarsening:" << std::endl;
  for(int i=ml-1; i>=0; i-=1)
  {
    error[i] = algorithm ( *grid , solution,-step, 1);
    if (i<ml-1) 
    {
      if( isLocallyAdaptive ) 
      {
        double eoc = log( error[i+1]/error[i]) / M_LN2; 
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

