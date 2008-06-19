#include <iostream>
#include <config.h>
#include <dune/common/stdstreams.cc>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

using namespace Dune;

#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/function/attachedfunction.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/space/dgspace/dgadaptmanager.hh>
#include <dune/fem/space/combinedspace/combinedadaptmanager.hh>

// #include "../leafindexset.hh"
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/space/lagrangespace.hh>

#if HAVE_GRAPE && GRIDDIM > 1 
#define USE_GRAPE 1
#else 
#define USE_GRAPE 0
#endif

#if USE_GRAPE && GRIDDIM > 1 
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/io/file/grapedataio.hh>

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
const int polOrd = POLORDER;

#ifndef GRIDDIM 
#define GRIDDIM dimworld 
#endif

#define GENERIC_ADAPT 1

//***********************************************************************
/*! L2 Projection of a function f: 
*/
//***********************************************************************

//! the index set we are using 
typedef DGAdaptiveLeafGridPart<GridType> GridPartType; 
//typedef HierarchicGridPart<GridType> GridPartType; 
//typedef AdaptiveLeafGridPart<GridType> GridPartType; 

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
// typedef MatrixFunctionSpace < double , double, GRIDDIM , 2,5 > FuncSpace;

//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh
typedef FunctionSpace < double , double, GRIDDIM , 5 > FuncSpace;
typedef FunctionSpace < double , double, GRIDDIM , 1 > SingleFuncSpace;
// typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
//    polOrd,CachingStorage> DiscreteFunctionSpaceType;
typedef DiscontinuousGalerkinSpace<SingleFuncSpace, GridPartType, 
  polOrd,CachingStorage> SingleDiscreteFunctionSpaceType;
typedef CombinedSpace<SingleDiscreteFunctionSpaceType,5,PointBased> 
   DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
//typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//typedef ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteFunctionSpaceType, DynamicVector< double > > > DiscreteFunctionType;
typedef AttachedDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

typedef AdaptationManager <GridType,
          RestrictProlongDefault<DiscreteFunctionType> > ADOperatorType;

// ***********************************************************
// the exact solution to the problem for EOC calculation 
class ExactSolution : public Function < FuncSpace , ExactSolution > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;
public:
  ExactSolution (const FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {}
 
  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 2.; // maximum of function is 2
    for(int i=0; i<DomainType::dimension; i++)
      ret *= sin(x[i]*(1.0 -x[i])*4.);
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};
 
// ********************************************************************
template <class DiscreteFunctionType, class FunctionType, int polOrd>
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

 public:
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
      FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;

    const FunctionSpaceType& space =  discFunc.space();

    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);

    Iterator it = space.begin();
    Iterator endit = space.end();

    // Get quadrature rule
    CachingQuadrature<GridPartType,0> quad(*it, 2*polOrd);

    for( ; it != endit ; ++it) {
      LocalFuncType lf = discFunc.localFunction(*it);
      const typename FunctionSpaceType::BaseFunctionSetType & baseset =
        lf.baseFunctionSet();
      const typename GridType::template Codim<0>::Entity::Geometry& 
                  itGeom = (*it).geometry();
      for(int qP = 0; qP < quad.nop(); qP++) 
      {
        f.evaluate(itGeom.global(quad.point(qP)), ret);
        for(int i=0; i<lf.numDofs(); i++) 
        {
          baseset.evaluate( i, quad[qP], phi );
          lf[i] += quad.weight(qP) * (ret * phi) ;
        }
      }
    }
  }
};
// calculates || u-u_h ||_L2
template <class DiscreteFunctionType>
class L2Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

public:
  template <int polOrd, class FunctionType>
  double norm (FunctionType &f, DiscreteFunctionType &discFunc,
      double time)
  {
    const typename DiscreteFunctionType::FunctionSpaceType
        & space = discFunc.space();

    typedef typename FunctionSpaceType::GridPartType GridPartType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typedef typename FunctionSpaceType::RangeType RangeType;

    RangeType ret (0.0);
    RangeType phi (0.0);

    double sum = 0.0;
    //LocalFuncType lf = discFunc.newLocalFunction();

    IteratorType it    = space.begin();
    IteratorType endit = space.end();

    // check whether grid is empty
    assert( it != endit );

    for(; it != endit ; ++it)
    {
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);
      LocalFuncType lf = discFunc.localFunction(*it);
      for(int qP = 0; qP < quad.nop(); qP++)
      {
        double det = (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
        lf.evaluate(quad[qP],phi);
        sum += det * quad.weight(qP) * SQR(ret[0] - phi[0]);
      }
    }
    return sqrt(sum);
  }
};
// ********************************************************************
void adapt(GridType& grid,
     DiscreteFunctionType& solution,int step) {
  typedef DiscreteFunctionType::FunctionSpaceType::Traits::IteratorType Iterator;
  const DiscreteFunctionType::FunctionSpaceType
    & space = solution.space();
  RestrictProlongDefault<DiscreteFunctionType> rp(solution);
  rp.setFatherChildWeight(DGFGridInfo<GridType>::refineWeight());

#if GENERIC_ADAPT 
  ADOperatorType adop(grid,rp);
#else 
  DofManagerType& dm = DofManagerFactoryType :: getDofManager( grid );  
#endif  

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

  
#if GENERIC_ADAPT
  for(int i=0; i<count; ++i)
  {
    Iterator it = space.begin();  
    Iterator endit = space.end();
    for(; it != endit ; ++it) {
      grid.mark(mark,it);
    } 
    adop.adapt();
    std::cout << message << std::endl;
  }
#else 
  for(int i=0; i<count; ++i)
  {
    Iterator it = space.begin();  
    Iterator endit = space.end();
    for(; it != endit ; ++it) {
      grid.mark(mark,it);
    } 
    grid.adapt(dm,rp);
    std::cout << message << "NOT GERNERIC!" << std::endl;
  }
#endif
  /*
  DofManagerType& dm = DofManagerFactoryType :: getDofManager( grid );  
  grid.preAdapt();
  dm.resizeForRestrict();
  grid.globalRefine(step);
  dm.resize();
  // dm.dofCompress();
  grid.postAdapt();
  */
}
// ********************************************************************
double algorithm (GridType& grid, DiscreteFunctionType& solution,
      int step,
      int turn )
{
  {
    const DiscreteFunctionSpaceType & space = solution.space();
    ExactSolution f ( space ); 
    L2Projection<DiscreteFunctionType, ExactSolution, polOrd>::
      project(f, solution);
    L2Error < DiscreteFunctionType > l2err;
    double new_error = l2err.norm<polOrd + 4> (f ,solution, 0.0);
    std::cout << "before ref." << new_error << "\n\n"; 
  }
  adapt(grid,solution,step);
#if USE_GRAPE
  // if Grape was found, then display last solution 
  if(0 && turn > 0) {
    std::cerr << "GRAPE 1" << std::endl;
    GrapeDataDisplay < GridType > grape(grid); 
    grape.dataDisplay( solution );
  }
#endif
  const DiscreteFunctionSpaceType & space = solution.space();
  ExactSolution f ( space ); 
  // calculation L2 error on refined grid
  // pol ord for calculation the error chould by higher than 
  // pol for evaluation the basefunctions 
  L2Error < DiscreteFunctionType > l2err;
  double error = l2err.norm<polOrd + 4> (f ,solution, 0.0);

#if USE_GRAPE
  // if Grape was found, then display last solution 
  if(0 && turn > 0) {
    std::cerr << "GRAPE 2" << std::endl;
    GrapeDataDisplay < GridType > grape(grid); 
    grape.dataDisplay( solution );
  }
#endif
  
  //! perform l2-projection to refined grid
  L2Projection<DiscreteFunctionType, ExactSolution, polOrd>::
    project(f, solution);
  double new_error = l2err.norm<polOrd + 4> (f ,solution, 0.0);
  std::cout << "\nL2 Error : " << error << " on new grid " << new_error << "\n\n";
#if USE_GRAPE
  // if Grape was found, then display last solution 
  if(0 && turn > 0) {
    std::cerr << "SIZE: " << solution.space().size() 
	      << " GRID: " << grid.size(0) << std::endl;
    std::cerr << "GRAPE 3" << std::endl;
    GrapeDataDisplay < GridType > grape(grid); 
    grape.dataDisplay( solution );
  }
#endif
  GrapeDataIO< GridType > dataio; 
  dataio.writeGrid( grid, xdr, "gridout", 0.0, turn );
  dataio.writeData( solution, xdr, "sol", turn );
  
  return error;
}


//**************************************************
//
//  main programm, run algorithm twice to calc EOC 
//
//**************************************************
int main (int argc, char **argv)
{
  if(argc != 2)
  {
    fprintf(stderr,"usage: %s <maxlevel> \n",argv[0]);
    exit(1);
  }
  int ml = atoi( argv[1] );
  std::vector<double> error(ml);

  char tmp[100]; 
  sprintf(tmp,"%ddgrid.dgf",GRIDDIM);
  GridPtr<GridType> gridptr(tmp);
  GridType* grid= gridptr.operator -> ();

  const int step = DGFGridInfo<GridType>::refineStepsForHalf();

  GridPartType part ( *grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  solution.clear();
  std::cout << "------------    Refining:" << std::endl;
  for(int i=0; i<ml; i+=1)
  {
    error[i] = algorithm ( *grid , solution, step, (i==ml-1));
    if (i>0) {
      double eoc = log( error[i-1]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  std::cout << "------------   Coarsening:" << std::endl;
  for(int i=ml-1; i>=0; i-=1)
  {
    error[i] = algorithm ( *grid , solution,-step, 1);
    if (i<ml-1) {
      double eoc = log( error[i+1]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  return 0;
}

