#include <config.h>

#include <iostream>
#include <dune/common/stdstreams.cc>

using namespace Dune;

#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/function/attachedfunction.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/dgspace/dgadaptmanager.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#if HAVE_GRAPE && GRIDDIM > 1 
#define USE_GRAPE 1
#else 
#define USE_GRAPE 0
#endif

#if USE_GRAPE && GRIDDIM > 1 
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/io/file/binarydataio.hh>
#include <dune/fem/io/parameter.hh>

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

//! define the type of discrete function we are using , see
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//typedef ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteFunctionSpaceType, DynamicVector< double > > > DiscreteFunctionType;
//typedef AttachedDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef DofManager< MyGridType > DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

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
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
      FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename GridPartType :: template Codim < 0 > :: 
      template Partition< All_Partition > :: IteratorType IteratorType;

    const FunctionSpaceType& space =  discFunc.space();

    if( polOrd < 0 )  polOrd = 2 * space.order() + 2 ;

    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);

    IteratorType it    = space.gridPart().template begin< 0, All_Partition > ();
    IteratorType endit = space.gridPart().template end< 0, All_Partition > ();

    // Get quadrature rule
    CachingQuadrature<GridPartType,0> quad(*it, polOrd);

    for( ; it != endit ; ++it) 
    {
      LocalFuncType lf = discFunc.localFunction(*it);
      const typename MyGridType::template Codim<0>::Entity::Geometry &itGeom
        = (*it).geometry();
      for( size_t qP = 0; qP < quad.nop(); ++qP )
      {
        f.evaluate(itGeom.global(quad.point(qP)), ret);
        ret *= quad.weight(qP) ;
        lf.axpy( quad[qP], ret );
      }
    }
  }
};

// calculates || u-u_h ||_L2
template <class DiscreteFunctionType>
class L2ErrorNoComm
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

public:
  template<class FunctionType >
  double norm (FunctionType &f, DiscreteFunctionType &discFunc, int polOrd = -1 )
  {
    const DiscreteFunctionSpaceType &space = discFunc.space();

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

    typedef typename GridPartType :: template Codim< 0 > :: 
      template Partition< All_Partition > :: IteratorType IteratorType ;

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    if( polOrd < 0 ) polOrd = 2*space.order() + 4 ; 

    RangeType ret (0.0);
    RangeType phi (0.0);

    double sum = 0.0;

    IteratorType it    = space.gridPart().template begin< 0, All_Partition > ();
    IteratorType endit = space.gridPart().template end< 0, All_Partition > ();

    for(; it != endit ; ++it)
    {
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);
      LocalFuncType lf = discFunc.localFunction(*it);
      for( size_t qP = 0; qP < quad.nop(); ++qP )
      {
        double det = (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)), ret);
        lf.evaluate(quad[qP],phi);
        RangeType diff = ret - phi ;
        sum += det * quad.weight(qP) * ( diff * diff );
      }
    }
    return sqrt(sum);
  }
};

void resetNonInterior( DiscreteFunctionType &solution )
{
  typedef DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
  typedef DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef DiscreteFunctionSpaceType :: GridPartType GridPartType ;
  typedef GridPartType :: Codim< 0 > :: Partition< All_Partition > :: IteratorType
    IteratorType;

  typedef IteratorType :: Entity  EntityType ;

  const DiscreteFunctionSpaceType& space = solution.space();

  IteratorType it    = space.gridPart().begin< 0, All_Partition > ();
  IteratorType endit = space.gridPart().end< 0, All_Partition > ();

  int count = 0;
  for( ; it != endit; ++ it ) 
  {
    const EntityType& entity = * it ;
    if( entity.partitionType() != InteriorEntity ) 
    {
      ++count ;
      solution.localFunction( entity ).clear();
    }
  }

  std::cout << "P[" << space.grid().comm().rank() << "]  reset " << count << " entities "  << std::endl;
            
}
 
// ********************************************************************
double algorithm ( MyGridType &grid, DiscreteFunctionType &solution, int step, int turn )
{
  ExactSolution f; 
  L2ErrorNoComm< DiscreteFunctionType > l2err;
  solution.clear();

  DGL2ProjectionAllPartitionNoComm :: project( f, solution );
  double new_error = l2err.norm(f ,solution);
  std::cout << "P[" << grid.comm().rank() << "]  start comm: " << new_error << std::endl; 

#if USE_GRAPE
  // if Grape was found, then display last solution 
  if(0 && turn > 0) {
    std::cerr << "GRAPE 1" << std::endl;
    GrapeDataDisplay < MyGridType > grape(grid); 
    grape.dataDisplay( solution );
  }
#endif

  // reset all non-interior entities, 
  // these should be restored during communication
  resetNonInterior( solution );

  // do communication 
  solution.communicate(); 

  // calculate l2 error again 
  double error = l2err.norm(f ,solution);
  std::cout << "P[" << grid.comm().rank() << "]  done comm: " << error << std::endl;

#if USE_GRAPE
  // if Grape was found, then display last solution 
  if(0 && turn > 0) {
    std::cerr << "GRAPE 2" << std::endl;
    GrapeDataDisplay< MyGridType > grape(grid); 
    grape.dataDisplay( solution );
  }
#endif

  if( std::abs( new_error - error ) > 1e-10 ) 
    DUNE_THROW(InvalidStateException,"Communication not working correctly");
  
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
  //ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  std::vector<double> error(ml);

  std::ostringstream gridFilenameStream;
  gridFilenameStream << GRIDDIM << "dgrid.dgf";
  GridPtr< MyGridType > gridptr ( gridFilenameStream.str() );
  MyGridType& grid = *gridptr ;

  const int step = DGFGridInfo< MyGridType >::refineStepsForHalf();

  GridPartType part ( grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  for(int i=0; i<ml; ++i )
  {
    if( grid.comm().rank() == 0) 
      std::cout << std::endl << "**** Start communication cycle " << i << "  ****" << std::endl;
    GlobalRefine :: apply( grid, step ); 
    error[ i ] = algorithm ( grid , solution, step, (i==ml-1));
  }
  return 0;
}
catch( const Dune :: Exception &exception )
{
  std :: cerr << exception << std :: endl;
  return 1;
}

