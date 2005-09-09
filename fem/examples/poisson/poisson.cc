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

#define SGRID 0

#if SGRID
#include <dune/grid/sgrid.hh>
#else
#include <dune/grid/albertagrid.hh>
#endif

#include <dune/fem/discreteoperatorimp.hh>
#include <dune/fem/lagrangebase.hh>
#include <dune/fem/dfadapt.hh>
//#include <dune/fem/discretefunction/adaptivefunction.hh>

#include "laplace.hh"

#include <dune/grid/common/leafindexset.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/grid/common/referenceelements.hh>

#include <dune/io/visual/grapedatadisplay.hh>

#include "../../solver/oemsolver/oemsolver.hh"

// laplace operator and L2 projection and error 
#include "laplace.cc"

using namespace Dune;

//***********************************************************************
/*! Poisson problem: 

  This is an example how to solve the equation on 
  \f[\Omega = (0,1)^2 \f]

  \f[ -\triangle u  = f \ \ \ in \Omega \f]
  \f[  \qquad u = 0  \ \ \ on  \partial\Omega \f]
  \f[ f(x,y) = 2 ( x + y - x^2 - y^2 ) \f]

  An exact solution to this problem is given by 
  \f[ u(x,y) = x ( 1 - x) y ( 1 - y ) \f]

  with the finite element method using lagrangian elements of polynor order 1.
*/
//***********************************************************************



#if SGRID
typedef SGrid     < dimp, dimw > GridType;
#else 
typedef AlbertaGrid< dimp, dimw > GridType;
#endif

// forward declaration 
class Tensor; 

/* define the types we are using */
//! dofmanger in not needed for DiscFuncArray, but we have to define it anyway 
// see dune/fem/dofmanager.hh
typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

typedef DefaultGridIndexSet<GridType,LevelIndex> IndexSetType;
//typedef DefaultGridIndexSet<GridType,GlobalIndex> IndexSetType;
//typedef DefaultGridPart<GridType,IndexSetType> GridPartType;
//typedef GridType :: LeafIndexSet IndexSetType;
//typedef LevelGridPart < GridType > GridPartType;
typedef DefaultGridPart<GridType,IndexSetType> GridPartType;

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

//! define the discrete laplace operator, see ./fem.cc
typedef LaplaceFEOp< DiscreteFunctionType, Tensor, 1 > LaplaceOperatorType;

//! define the inverse operator we are using to solve the system 
// see dune/fem/inverseoperators.hh 
//typedef CGInverseOp < DiscreteFunctionType, LaplaceOperatorType >    InverseOperatorType;
/****************************************/
// or ../../solvers/oemsolver/oemsolvers.hh
//typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;

//! define the type of mapping which is used by inverseOperator 
typedef Mapping<double ,double,DiscreteFunctionType,DiscreteFunctionType > MappingType;

// right hand side of governing problem 
class RHSFunc : public Function < FuncSpace , RHSFunc > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::DomainType DomainType;
 
public:
  RHSFunc (FuncSpace &f)
    : Function <  FuncSpace , RHSFunc > (f) {}
   
  // f(x,y) = 2 * (x + y - x*x -y*y )
  void evaluate (const DomainType & x , RangeType & ret) 
  {
    ret = 0.0;
    for(int i=0; i<DomainType::dimension; i++)
      ret += (x[i] -  SQR(x[i])); 
    ret *= 2.0;
  }
};

//! the exact solution to the problem for EOC calculation 
class ExactSolution : public Function < FuncSpace , ExactSolution > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;
public:
  ExactSolution (FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {}
 
  //! u(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret) 
  {
    ret = 1.0;
    for(int i=0; i<DomainType::dimension; i++)
      ret *= x[i]*(1.0 -x[i]);
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) 
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
  typedef typename EntityType::IntersectionIterator NeighIt;
  typedef typename NeighIt::BoundaryEntity BoundaryEntityType;

  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::GridType GridType;

  const FunctionSpaceType & space = rhs.getFunctionSpace();
    
  typedef typename DiscreteFunctionType::DofIteratorType DofIterator;
  DofIterator dit = rhs.dbegin();
      
  NeighIt endit = en.iend();
  for(NeighIt it = en.ibegin(); it != endit; ++it)
  {
    if(it.boundary())
    {
      typedef typename EntityType :: ctype coordType; 
      enum { dim = EntityType :: dimension };
      GeometryType t = en.geometry().type();
      
      if( (t == simplex) || (t == triangle) || (t == tetrahedron ) )
      {
        static ReferenceSimplex< coordType, dim > refElem; 
        int face = it.numberInSelf();
        int novx = refElem.size( face, 1 , dim );
        assert( novx == dim );
        for(int j=0; j< novx ; j++)
        {
          int vx = refElem.subentity(face,1, j , dim );
          int row = space.mapToGlobal( en , vx );
          dit[row] = 0.0;
        }
      }
      if( en.geometry().type() == cube )
      {
        static ReferenceCube < coordType, dim > refElem; 
        int face = it.numberInSelf();
        int novx = refElem.size( face, 1 , dim );
        for(int j=0; j< novx ; j++)
        {
          int vx = refElem.subentity(face,1, j , dim );
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

   IndexSetType iset ( grid , grid.maxlevel () );
   //const IndexSetType & iset = grid.leafIndexSet();
   std::cout << "\nSolving for " << iset.size(dimp) << " number of unkowns. \n\n";
   
   GridPartType part ( grid, iset );

   FuncSpaceType linFuncSpace ( part );
   DiscreteFunctionType solution ( "sol", linFuncSpace );
   solution.clear();
   DiscreteFunctionType rhs ( "rhs", linFuncSpace );
   rhs.clear();
      
   RHSFunc f ( linFuncSpace ); 
    
   LaplaceOperatorType laplace ( linFuncSpace , LaplaceOperatorType::ASSEMBLED);
   
   //! build right hand side, does not allocate b!
   L2Projection < DiscreteFunctionType > l2pro;
   const int polOrd = 2;
   l2pro.project<polOrd> ( f , rhs );
    
   { 
     typedef FuncSpaceType :: IteratorType IteratorType; 
     // set Dirichlet Boundary to zero 
     IteratorType endit  = linFuncSpace.end();
     for(IteratorType it = linFuncSpace.begin(); it != endit; ++it ) 
     {
       boundaryTreatment ( *it , rhs );
     }
   }

   //laplace.print();
   //rhs.print(std::cout);
    
   bool verbose = false; 
   InverseOperatorType cg ( laplace, 1E-3 , 1E-6 , 20000 , verbose );
     
   // solve linear system with cg 
   cg(rhs,solution);

   // calculation L2 error 
   ExactSolution u ( linFuncSpace ); 
   L2Error < DiscreteFunctionType > l2err;

   // pol ord for calculation the error chould by higher than 
   // pol for evaluation the basefunctions 
   double error = l2err.norm<polOrd + 2> (u ,solution, 0.0);
   std::cout << "\nL2 Error : " << error << "\n\n";

   
   if(turn > 0)
   {
     //GrapeDataDisplay < GridType , DiscreteFunctionType > grape(grid); 
     //grape.dataDisplay( solution );
   }
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
  //std::string macroGridName ( "2drefel.al" );
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

