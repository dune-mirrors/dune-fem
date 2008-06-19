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

#define SGRID 1
#define AGRID 0

using namespace Dune;

#if SGRID
#include <dune/grid/sgrid.hh>
typedef SGrid  < dimp, dimw > GridType;
#endif

#if AGRID  
#include <dune/grid/albertagrid.hh>
typedef AlbertaGrid< dimp, dimw > GridType;
#endif

//- Dune includes 
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/grid/common/referenceelements.hh>

#include "../../operator/discreteoperatorimp.hh"
#include "../../space/lagrangespace.hh"
#include "../../discretefunction/dfadapt.hh"

#include <dune/fem/space/common/adaptiveleafgridpart.hh>
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include "../../operator/inverseoperators.hh"
#include "../../solver/oemsolver/oemsolver.hh"

#include <dune/fem/misc/l2error.hh>

// massmatrix and L2 projection and error 
#include "../poisson/laplace.cc"
#include "massmatrix.hh"

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
const int polOrd = 2;

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

// forward declaration 
class Tensor; 

//! the index set we are using 
//typedef LevelGridPart < GridType > GridPartType;
typedef LeafGridPart<GridType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
typedef FunctionSpace < double , double, dimp , 1 > FuncSpace;

//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace < FuncSpace , GridPartType , 1 > FuncSpaceType ;

//! define the type of discrete function we are using , see
//! dune/fem/discfuncarray.hh
typedef DFAdapt < FuncSpaceType > DiscreteFunctionType;
//typedef AdaptiveDiscreteFunction < FuncSpaceType > DiscreteFunctionType;
//typedef DiscFuncArray < FuncSpaceType > DiscreteFunctionType;

//! define the discrete massmatrix operator, see ./fem.cc
//typedef LaplaceFEOp< DiscreteFunctionType, Tensor, 1 > LaplaceOperatorType;
//typedef LaplaceFEOp< DiscreteFunctionType, Tensor, 1 > LaplaceOperatorType;
typedef MassMatrixFEOp<DiscreteFunctionType,Tensor, polOrd + 2>  LaplaceOperatorType;

//! define the inverse operator we are using to solve the system 
// see dune/fem/inverseoperators.hh 
typedef CGInverseOp < DiscreteFunctionType, LaplaceOperatorType >    InverseOperatorType;
/****************************************/
// or ../../solvers/oemsolver/oemsolvers.hh
//typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;

//! the exact solution to the problem for EOC calculation 
class ExactSolution : public Function < FuncSpace , ExactSolution > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;
public:
  ExactSolution (FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {}
 
  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret) const
  {
    ret = 1.0;
    for(int i=0; i<DomainType::dimension; i++)
      ret *= x[i]*(1.0 -x[i]);
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};
 
// diffusion coefficient for this problem the id 
class Tensor : public Function < FunctionSpace < double , double, dimp , 1 >, Tensor >
{
  typedef FunctionSpace< double , double, dimp , 1 >  FuncSpace;
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::DomainType DomainType;

public:
  // Constructor 
  Tensor (FuncSpace &f)
    : Function < FuncSpace , Tensor > ( f )  { } ;
  // eval Tensor 
  void evaluate (int i, int j, const DomainType & x1 , RangeType & ret)
  {
    evaluate(x1,0.0,ret);
  }
  void evaluate (const DomainType & x1 , RangeType & ret)
  {
    evaluate(x1,0.0,ret);
  }
  void evaluate (const DomainType & x1 ,double time, RangeType & ret)
  {
    ret[0] = 1.0;
  }
};//end class Tensor


//! set the dirichlet points to zero 
template <class EntityType, class DiscreteFunctionType> 
void boundaryTreatment ( const EntityType & en ,  DiscreteFunctionType &rhs )
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::GridPartType GridPartType;
  typedef typename FunctionSpaceType::GridType GridType;

  typedef typename GridPartType :: IntersectionIteratorType NeighIt;

  const FunctionSpaceType & space = rhs.space();
  const GridPartType& gridPart = space.gridPart();
    
  typedef typename DiscreteFunctionType::DofIteratorType DofIterator;
  DofIterator dit = rhs.dbegin();
      
  NeighIt endit = gridPart.iend(en);
  for(NeighIt it = gridPart.ibegin(en); it != endit; ++it)
  {
    if(it.boundary())
    {
      typedef typename EntityType :: ctype coordType; 
      enum { dim = EntityType :: dimension };
      GeometryType t = en.geometry().type();
      
      if( t.isSimplex() )
      {
        static ReferenceSimplex< coordType, dim > refElem; 
        int face = it.numberInSelf();
        int novx = refElem.size( face, 1 , dim );
        assert( novx == dim );
        for(int j=0; j< novx ; j++)
        {
          int vx = refElem.subEntity(face,1, j , dim );
          int row = space.mapToGlobal( en , vx );
          dit[row] = 0.0;
        }
      }
      if( t.isCube() )
      {
        static ReferenceCube < coordType, dim > refElem; 
        int face = it.numberInSelf();
        int novx = refElem.size( face, 1 , dim );
        for(int j=0; j< novx ; j++)
        {
          int vx = refElem.subEntity(face,1, j , dim );
          int row = space.mapToGlobal( en , vx);
          dit[row] = 0.0;
        }
      }
    }
  }
}

double algorithm (const char * filename , int maxlevel, int turn )
{
  // we dont not use all levels of grid for calculation, only maxlevel
#if SGRID
   // this leads to the same number of points for SGrid and AlbertGrid
   int n[dimp];
   double h[dimp];
   for(int i=0; i<dimp; i++)  { n[i] = 2; h[i] = 1.0; }

   GridType grid ((int *) &n, (double *) &h );
#else
   GridType grid ( filename );
#endif

   grid.globalRefine (maxlevel);

   GridPartType part ( grid );

   FuncSpaceType linFuncSpace ( part );
   std::cout << "\nSolving for " << linFuncSpace.size() << " number of unkowns. \n\n";
   
   DiscreteFunctionType solution ( "sol", linFuncSpace );
   solution.clear();
   DiscreteFunctionType rhs ( "rhs", linFuncSpace );
   rhs.clear();
      
   ExactSolution f ( linFuncSpace ); 
   L2Error < DiscreteFunctionType > l2err;
    
   LaplaceOperatorType massmatrix ( linFuncSpace , LaplaceOperatorType::ASSEMBLED);
   
   //! build right hand side, does not allocate b!
   L2Projection < DiscreteFunctionType > l2pro;
   l2pro.project<polOrd + 2> ( f , rhs );
   
   { 
     typedef FuncSpaceType :: IteratorType IteratorType; 
     // set Dirichlet Boundary to zero 
     IteratorType endit  = linFuncSpace.end();
     for(IteratorType it = linFuncSpace.begin(); it != endit; ++it ) 
     {
       boundaryTreatment ( *it , rhs );
     }
   }

   //massmatrix.print();
   //rhs.print(std::cout);
    
   bool verbose = false; 
   InverseOperatorType cg ( massmatrix, 1E-3 , 1E-6 , 20000 , verbose );
     
   // solve linear system with cg 
   cg(rhs,solution);

   // calculation L2 error 
   // pol ord for calculation the error chould by higher than 
   // pol for evaluation the basefunctions 
   double error = l2err.norm( f ,solution);
   std::cout << "\nL2 Error : " << error << "\n\n";
  
#if HAVE_GRAPE
   // if Grape was found, then display last solution 
   if(turn > 0)
   {
     GrapeDataDisplay < GridType > grape(grid); 
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
int main (int argc, char **argv)
{
  if(argc != 2)
  {
    fprintf(stderr,"usage: %s <maxlevel> \n",argv[0]);
    exit(1);
  }
  int ml = atoi( argv[1] );
  double error[2];
  char tmp[16]; sprintf(tmp,"%d",dimp);
  std::string macroGridName (tmp); 
  macroGridName += "dgrid.al";

#if SGRID 
  const int step = 1;
#else 
  const int step = 2;
#endif
  
  ml -= step;
  if(ml < 0) ml = 0;
  for(int i=0; i<2; i++)
  {
    error[i] = algorithm ( macroGridName.c_str() ,  ml , i);
    ml += step ;
  }
  double eoc = log( error[0]/error[1]) / M_LN2; 
  std::cout << "EOC = " << eoc << " \n";
  return 0;
}

