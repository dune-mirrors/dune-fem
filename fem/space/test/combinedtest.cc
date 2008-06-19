#include <iostream>
#include <config.h>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
static const int dimw = dimworld;

#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/space/combinedspace/combinedspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/combinedfunction.hh>
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
// double -> Double

// the exact solution to the problem for EOC calculation 
template <class FuncSpace>
class ExactSolution : public Function < FuncSpace , ExactSolution<FuncSpace> > 
{
  typedef Function < FuncSpace , ExactSolution<FuncSpace> > BaseType;
  typedef typename FuncSpace::RangeType RangeType;
  typedef typename FuncSpace::RangeFieldType RangeFieldType;
  typedef typename FuncSpace::DomainType DomainType;
public:
  ExactSolution (const FuncSpace &f,int shift=0) :
    BaseType ( f ), shift_(shift) {}
 
  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 2.; // maximum of function is 2
    for (int j=0;j<RangeType::dimension; j++) 
      for(int i=0; i<DomainType::dimension; i++)
	ret[j] += pow(x[i]*(1.0 -x[i])*4.,double(j+1+shift_));
  }
  void evaluate (const DomainType & x , RangeFieldType time , 
		 RangeType & ret) const
  {
    evaluate ( x , ret );
  }
  int shift_;
};
 
// ********************************************************************
template <class DiscreteFunctionType>
class L2Projection
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType 
          DiscreteFunctionSpaceType;

 public:
  template <class FunctionType>
  static void project (const FunctionType &f, 
         DiscreteFunctionType &discFunc, int polOrd) 
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
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType 
	BaseFunctionSetType ; 
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
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType 
          DiscreteFunctionSpaceType;
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
// ********************************************************************
// ********************************************************************
const int RANGE = 3;
typedef FunctionSpace < GridType :: ctype, double , 
			dimw , RANGE > FuncSpace;
typedef FunctionSpace < GridType :: ctype, double , 
			dimw , 1 > SingleFuncSpace;
typedef DiscontinuousGalerkinSpace<
  SingleFuncSpace,GridPartType, 
  polOrd,CachingStorage> SingleDiscreteFunctionSpaceType;
typedef AdaptiveDiscreteFunction < 
  SingleDiscreteFunctionSpaceType > SingleDiscreteFunctionType;

//! Get the Dofmanager type
typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> 
        DofManagerFactoryType;

template <class DiscreteFunctionType>
double algorithm (GridType& grid, 
		  DiscreteFunctionType& solution  , int turn )
{
  typedef typename 
    DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
   //! perform l2-projection
   typedef typename DiscreteFunctionSpaceType :: 
     ContainedSpaceType SubDFSType;
   typedef typename DiscreteFunctionType :: 
     SubDiscreteFunctionType SubDFType;
   typedef typename DiscreteFunctionSpaceType :: 
     RangeType RangeType; 
   typedef typename SubDFSType :: RangeType SubRangeType;
   RangeType error(0);
  const DiscreteFunctionSpaceType& 
    linFuncSpace = solution.space();
  ExactSolution<DiscreteFunctionSpaceType> f (linFuncSpace); 
   L2Projection<DiscreteFunctionType>:: project(f, solution);
   for (int method=0;method<2;method++) {
     // calculation L2 error 
     // pol ord for calculation the error chould by higher than 
     // pol for evaluation the basefunctions 
     L2Error < DiscreteFunctionType > l2err;
     error = l2err.norm(f ,solution, 0.0);
     for(int i=0; i<RangeType::dimension; ++i)
       std::cout << "\nL2 Error["<<i<<"] : " << error[i] << "\n\n";
     for (int i=0;i<RangeType::dimension; ++i) {
       SubDFType& sol0 = solution.subFunction(i);
       ExactSolution<SubDFSType> f0 (linFuncSpace.containedSpace(),i ); 
       L2Error < SubDFType > l2err0;
       SubRangeType error0 = l2err0.norm(f0,sol0,0.0);
       std::cout << "\n L2 SubError[" << i << "]-Error[" << i << "]  : " 
                 << error0[0]-error[i] << "\n\n";
     }
     solution.clear();
     for (int i=0;i<RangeType::dimension; i+=1) {
       SubDFType& sol0 = solution.subFunction(i);
       ExactSolution<SubDFSType> f0 (sol0.space(),i ); 
       L2Projection<SubDFType>:: project(f0, sol0);
     }
   }
   
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

  std :: cout << "Polynomial Order: " << polOrd << std :: endl;

  typedef CombinedDiscreteFunction<
    SingleDiscreteFunctionType,RANGE > DiscreteFunctionType1;
  DiscreteFunctionType1* solution1;
  {
    GridPtr<GridType> gridptr(macroGridName);
    GridType& grid=*gridptr;
    const int step = Dune::DGFGridInfo<GridType>::
      refineStepsForHalf();
    GridPartType part ( grid );
    SingleDiscreteFunctionSpaceType singFuncSpace ( part );
    SingleDiscreteFunctionType singSol("sol",singFuncSpace);
    solution1 = new DiscreteFunctionType1( singSol );
    solution1->clear();
    for(int i=0; i<ml; i+=step) {
      grid.globalRefine(step);
      DofManagerType& dm = DofManagerFactoryType :: 
	      getDofManager( grid );
      dm.resize();
      error[i] = algorithm ( grid , *solution1 , i==ml-1);
      if (i>0) {
	double eoc = log( error[i-step]/error[i]) / M_LN2; 
	std::cout << "EOC = " << eoc << " \n";
      }
    }
  }
  typedef CombinedSpace
     <SingleDiscreteFunctionSpaceType,RANGE,VariableBased> 
      DiscreteFunctionSpaceType2;
  typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType2>
      DiscreteFunctionType2;
  DiscreteFunctionType2* solution2;
  {
    GridPtr<GridType> gridptr(macroGridName);
    GridType& grid=*gridptr;
    const int step = Dune::DGFGridInfo<GridType>::
      refineStepsForHalf();
    GridPartType part ( grid );
    DiscreteFunctionSpaceType2 funcSpace(part);
    solution2 = new DiscreteFunctionType2 ( "sol",funcSpace );
    solution2->clear();
    for(int i=0; i<ml; i+=step) {
      grid.globalRefine(step);
      DofManagerType& dm = DofManagerFactoryType :: 
	getDofManager( grid );
      dm.resize();
      error[i] = algorithm ( grid , *solution2 , i==ml-1);
      if (i>0) {
	double eoc = log( error[i-step]/error[i]) / M_LN2; 
	std::cout << "EOC = " << eoc << " \n";
      }
    }
  }
  (*solution1) -= (*solution2);
  std::cout << "Difference of function: " 
	    << solution1->scalarProductDofs(*solution1) 
            << std::endl;
  delete solution1;
  delete solution2;
  delete [] error;
  return 0;
}

