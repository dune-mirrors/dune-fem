#include <iostream>
#include <config.h>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
static const int dimw = dimworld;

#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/fem/space/common/adaptiveleafgridpart.hh> 
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/misc/double.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif
using namespace Dune;

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
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

//! the index set we are using 
typedef HierarchicGridPart<GridType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
//typedef MatrixFunctionSpace < double , double, dimw , 3,5 > FuncSpace;
typedef FunctionSpace < GridType :: ctype, Double , dimw , 2 > FuncSpace;

//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
	polOrd,CachingStorage> DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
//! dune/fem/discfuncarray.hh
typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

//! Get the Dofmanager type
typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;


// the exact solution to the problem for EOC calculation 
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
    for (int j=0;j<RangeType::dimension; j++) 
      for(int i=0; i<DomainType::dimension; i++)
	ret[j] *= pow(x[i]*(1.0 -x[i])*4.,double(j+1));
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};
 
// ********************************************************************
template <class DiscreteFunctionType>
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;

 public:
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc, int polOrd) 
  {
    typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;

    const DiscreteFunctionSpaceType& space =  discFunc.space();

    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typename DiscreteFunctionSpaceType::RangeType ret (0.0);
    typename DiscreteFunctionSpaceType::RangeType phi (0.0);

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      // Get quadrature rule
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);

      LocalFuncType lf = discFunc.localFunction(*it);

      //! Note: BaseFunctions must be ortho-normal!!!!
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
      const BaseFunctionSetType & baseset =
        lf.baseFunctionSet();

      const typename GridType::template Codim<0>::Entity::Geometry& 
        itGeom = (*it).geometry();
     
      const int quadNop = quad.nop();
      const int numDofs = lf.numDofs();
      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        f.evaluate(itGeom.global(quad.point(qP)), ret);
        for(int i=0; i<numDofs; ++i) {
          baseset.evaluate( i, quad[qP], phi );
          lf[i] += quad.weight(qP) * (ret * phi) ;
        }
      }
    }
  }
  
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) 
  {
    const DiscreteFunctionSpaceType& space =  discFunc.space();
    int polOrd = 2 * space.order();
    project(f,discFunc,polOrd);
  }
};


// calculates || u-u_h ||_L2
template <class DiscreteFunctionType>
class L2Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

public:
  template <class FunctionType>
  RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
      double time, int polOrd) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.space();

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    RangeType ret (0.0);
    RangeType phi (0.0);

    RangeType error(0.0);

    enum { dimRange = DiscreteFunctionSpaceType :: DimRange };

    IteratorType endit = space.end();
    for(IteratorType it = space.begin(); it != endit ; ++it)
    {
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);
      LocalFuncType lf = discFunc.localFunction(*it);
      const int quadNop = quad.nop();
      for(int qP = 0; qP < quadNop; ++qP)
      {
        double weight = quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
        lf.evaluate(quad[qP],phi);

        for(int i=0; i< dimRange; ++i)
          error[i] += weight * SQR(ret[i] - phi[i]);
      }
    }
    
    for(int i=0; i< dimRange; ++i) 
    {
      error[i] = sqrt(error[i]);
    }
    
    return error;
  }

  template <class FunctionType>
  RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
      double time) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.space();
    int polOrd = 2 * space.order() + 2;
    return norm(f,discFunc,time,polOrd);
  }
};
// ********************************************************************
double algorithm (GridType& grid, DiscreteFunctionType& solution  , int turn )
{
   GridPartType part ( grid );
   DiscreteFunctionSpaceType linFuncSpace ( part );
   ExactSolution f ( linFuncSpace ); 
   L2Error < DiscreteFunctionType > l2err;
       
   //! perform l2-projection
   L2Projection<DiscreteFunctionType>::
     project(f, solution);

   // calculation L2 error 
   // pol ord for calculation the error chould by higher than 
   // pol for evaluation the basefunctions 
   typedef DiscreteFunctionSpaceType :: RangeType RangeType; 
   RangeType error = l2err.norm(f ,solution, 0.0);
   for(int i=0; i<RangeType::dimension; ++i)
     std::cout << "\nL2 Error["<<i<<"] : " << error[i] << "\n\n";

#if HAVE_GRAPE
   // if Grape was found, then display last solution 
   if(0 && turn > 0)
   {
     GrapeDataDisplay < GridType > grape(part); 
     grape.dataDisplay( solution );
   }
#endif
   
   return sqrt(error*error);
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
  char tmp[16]; sprintf(tmp,"%d",dimw);
  std::string macroGridName (tmp); 
  macroGridName += "dgrid.dgf";

  GridPtr<GridType> gridptr(macroGridName);
  GridType& grid=*gridptr;
  const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();

  GridPartType part ( grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  solution.clear();
  
  for(int i=0; i<ml; i+=step)
  {
    grid.globalRefine(step);
    DofManagerType& dm = DofManagerFactoryType :: getDofManager( grid );
    dm.resize();
    error[i] = algorithm ( grid , solution , i==ml-1);
    if (i>0) 
    {
      double eoc = log( error[i-step]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  delete [] error;
  return 0;
}

