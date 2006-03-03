#ifndef DIM_OF_WORLD 
static const int dimw = 2;
#else
static const int dimw = DIM_OF_WORLD;
#endif

#ifndef DIM
static const int dimp = 2;
#else
static const int dimp = DIM;
#endif

#include <iostream>
#include <config.h>
#include <dune/common/stdstreams.cc>

#define AGRID 1

using namespace Dune;

#if AGRID  
#include <dune/grid/albertagrid.hh>
typedef AlbertaGrid< dimp, dimw > GridType;
#endif

#include <dune/fem/discreteoperatorimp.hh>
#include <dune/fem/lagrangebase.hh>
#include <dune/fem/dfadapt.hh>
#include "../../space/dgspace.hh"
#include "../../quadrature/cachequad.hh"
#include "../../space/adaptoperator.hh"

#include <dune/grid/common/leafindexset.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/grid/common/referenceelements.hh>

#if HAVE_GRAPE
#include <dune/io/visual/grapedatadisplay.hh>
#endif

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
const int polOrd = POLORDER;

//***********************************************************************
/*! L2 Projection of a function f: 
*/
//***********************************************************************

//! the index set we are using 
typedef DefaultGridIndexSet<GridType,GlobalIndex> IndexSetType;
// typedef AdaptiveLeafIndexSet<GridType> IndexSetType;
typedef DefaultGridPart<GridType,IndexSetType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
typedef FunctionSpace < double , double, dimp , 1 > FuncSpace;

//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
	polOrd,CachingStorage> DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
//! dune/fem/discfuncarray.hh
typedef DFAdapt < DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

typedef AdaptOperator<GridType,
		      RestProlOperator<DiscreteFunctionType> > ADOperatorType;

// ***********************************************************
//! the exact solution to the problem for EOC calculation 
class ExactSolution : public Function < FuncSpace , ExactSolution > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;
public:
  ExactSolution (FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {}
 
  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 2.; // maximum of function is 2
    for(int i=0; i<DomainType::dimension; i++)
      ret *= x[i]*(1.0 -x[i])*4.;
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
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;

    const FunctionSpaceType& space =  discFunc.getFunctionSpace();

    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);

    Iterator it = space.begin();
    Iterator endit = space.end();

    // Get quadrature rule
    CachingQuadrature<GridType,0> quad(*it, 2*polOrd+1);

    for( ; it != endit ; ++it) {
      LocalFuncType lf = discFunc.localFunction(*it);
      const typename FunctionSpaceType::BaseFunctionSetType & baseset =
        lf.getBaseFunctionSet();
      for(int i=0; i<lf.numDofs(); i++) {
        for(int qP = 0; qP < quad.nop(); qP++) {
          f.evaluate((*it).geometry().global(quad.point(qP)), ret);
          baseset.eval(i,quad,qP,phi);
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
        & space = discFunc.getFunctionSpace();

    typedef typename FunctionSpaceType::GridType GridType;
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
      CachingQuadrature<GridType,0> quad(*it, polOrd);
      LocalFuncType lf = discFunc.localFunction(*it);
      for(int qP = 0; qP < quad.nop(); qP++)
      {
        double det = (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
        lf.evaluate((*it),quad,qP,phi);
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
    & space = solution.getFunctionSpace();
  Iterator it = space.begin();  

  RestProlOperator<DiscreteFunctionType> rp(solution,it->geometry().type());
  ADOperatorType adop(grid,rp);

  Iterator endit = space.end();

  for(; it != endit ; ++it) {
    grid.mark(step,it);
  }
  adop.adapt();
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
  adapt(grid,solution,step);
  IndexSetType iset ( grid );
  GridPartType part ( grid, iset );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  ExactSolution f ( linFuncSpace ); 
  // calculation L2 error on refined grid
  // pol ord for calculation the error chould by higher than 
  // pol for evaluation the basefunctions 
  L2Error < DiscreteFunctionType > l2err;
  double error = l2err.norm<polOrd + 4> (f ,solution, 0.0);

#if HAVE_GRAPE
  // if Grape was found, then display last solution 
  if(turn > 0) {
    GrapeDataDisplay < GridType > grape(grid); 
    grape.dataDisplay( solution );
  }
#endif
  
  //! perform l2-projection to refined grid
  L2Projection<DiscreteFunctionType, ExactSolution, polOrd>::
    project(f, solution);
  double new_error = l2err.norm<polOrd + 4> (f ,solution, 0.0);
  std::cout << "\nL2 Error : " << error << " on new grid " << new_error << "\n\n";
  
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
  double* error = new double[ml];
  char tmp[16]; sprintf(tmp,"%d",dimp);
  std::string macroGridName (tmp); 
  macroGridName += "dgrid.al";

#if AGRID 
  const int step = 1;
#endif
#if AGRID
  GridType grid ( macroGridName.c_str() );
#endif

  IndexSetType iset ( grid );
  GridPartType part ( grid, iset );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  solution.clear();
  std::cout << "------------    Refining:" << std::endl;
  for(int i=0; i<ml; i+=step)
  {
    error[i] = algorithm ( grid , solution, step, (i==ml-1));
    if (i>0) {
      double eoc = log( error[i-step]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  /*
  std::cout << "------------   Coarsening:" << std::endl;
  for(int i=ml-1; i>=0; i-=step)
  {
    error[i] = algorithm ( grid , solution,-step, 1);
    if (i<ml-1) {
      double eoc = log( error[i+step]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  */
  delete [] error;
  return 0;
}

