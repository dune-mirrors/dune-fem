/**************************************************************************
**       Title: poisson.cc
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL Author: robertk 
** Description: File demonstrating a simple numerics problem on arbitrary
**              grids: poisson-problem with known solution is given
**              and compared with numerical solution (EOC)
**              Dune grid parser is used.
**
**              For changing grid-types, compile with
**
**              make clean
**
**              and one of the following
**           
**              make GRIDTYPE=YASPGRID       (default)
**                    -> compiles and works correctly
**              make
**                    -> compiles and works correctly
**              make GRIDTYPE=ALBERTAGRID
**                    -> compiles and works correctly
**              make GRIDTYPE=SGRID
**                    -> compiles and works correctly
**              make GRIDTYPE=ALUGRID_SIMPLEX
**                    -> compiles and runs, but EOCs are too low
** 
**
**-------------------------------------------------------------------------
**
**  $Log$
**  Revision 1.17  2006/09/25 12:30:10  haasdonk
**  now working for all grids except UG
**
**  Revision 1.16  2006/09/25 11:39:15  haasdonk
**  adjusting to dune-grid-parser
**
**
**************************************************************************/

#include <iostream>
#include <config.h>
#include <dune/common/stdstreams.cc>
#include <dune/grid/io/file/dgfparser/gridtype.hh>

using namespace Dune;

// #if SGRID
// #include <dune/grid/sgrid.hh>
// static const int dimworld = GRIDDIM;
// static const int dimworld = dimworld;
// typedef SGrid  < dimworld, dimworld > GridType;
// static const int refStepsForHalf = 1;
// #endif

// #if AGRID  
// #include <dune/grid/albertagrid.hh>
// static const int dimworld = DUNE_WORLD_DIM;
// static const int dimworld = DUNE_PROBLEM_DIM;

// typedef AlbertaGrid< dimworld, dimworld > GridType;
// static const int refStepsForHalf = dimworld;
// #endif

// #if BGRID  

// //#include <dune/grid/alu3dgrid/includecc.cc>
// //#include <dune/grid/alu3dgrid.hh>
// //#include <dune/grid/alugrid/3d/grid.hh>
// //static const int dimworld = 3;
// //static const int dimworld = 3;
// //typedef ALU3dGrid < dimworld, dimworld , tetra > GridType;

// #include <dune/grid/alugrid.hh>
// static const int dimworld = DUNE_PROBLEM_DIM;
// static const int dimworld = DUNE_PROBLEM_DIM;
// typedef ALUSimplexGrid < dimworld, dimworld > GridType;
// static const int refStepsForHalf = 1;
// #endif

//- Dune includes 
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/referenceelements.hh>

//- local inlcudes 
#include "../../operator/discreteoperatorimp.hh"
#include "../../space/lagrangespace.hh"
#include "../../discretefunction/dfadapt.hh"
#include "../../discretefunction/adaptivefunction.hh"

#include "laplace.hh"

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include "../../operator/inverseoperators.hh"
#include "../../solver/oemsolver/oemsolver.hh"

// laplace operator and L2 projection and error 
#include "laplace.cc"

//***********************************************************************
/*! Poisson problem: 

  This is an example how to solve the equation on 
  \f[\Omega = (0,1)^dimworld \f]

  \f[ -\triangle u  = f \ \ \ in \Omega \f]
  \f[  \qquad u = 0  \ \ \ on  \partial\Omega \f]
  \f[ f(x,y,z) = 2 (x-x^2) (y-y^2) +
                 2 (z-z^2) (y-y^2) +    
                 2 (x-x^2) (z-z^2)
  \f]

  An exact solution to this problem is given by 
  \f[ u(x,y,z) = ( x - x^2 ) ( y - y^2 ) ( z - z^2 ) \f]

  with the finite element method using lagrangian elements of polynor order 1.
*/
//***********************************************************************

// forward declaration 
class Tensor; 

//! the index set we are using 
//typedef LevelGridPart < GridType > GridPartType;
typedef LeafGridPart<GridType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
typedef FunctionSpace < double , double, dimworld , 1 > FuncSpace;

//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace < FuncSpace , GridPartType , 1 > FuncSpaceType ;

//! define the type of discrete function we are using , see
//! dune/fem/discfuncarray.hh
//typedef DFAdapt < FuncSpaceType > DiscreteFunctionType;
typedef AdaptiveDiscreteFunction < FuncSpaceType > DiscreteFunctionType;

//! define the discrete laplace operator, see ./fem.cc
typedef LaplaceFEOp< DiscreteFunctionType, Tensor, 1 > LaplaceOperatorType;

//! define the inverse operator we are using to solve the system 
// see dune/fem/inverseoperators.hh 
//typedef CGInverseOp < DiscreteFunctionType, LaplaceOperatorType >    InverseOperatorType;
/****************************************/
// or ../../solvers/oemsolver/oemsolvers.hh
typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
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
   
  //  f(x,y,z) = 2 (x-x^2) (y-y^2) +
  //             2 (z-z^2) (y-y^2) +              
  //             2 (x-x^2) (z-z^2)
  void evaluate (const DomainType & x , RangeType & ret) 
  {
    enum { dim = DomainType::dimension };
    ret = 0.0;
    for(int i=0; i<dim; i++)
    { 
      RangeType tmp = 2.0;
      for(int j=1; j<dim; j++)
      {
        int idx = (i+j) % dim;
        tmp *= (x[idx] - SQR(x[idx]));
      }
      ret += tmp;
    }
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
 
  //! u(x,y,z) = (x-x^2)*(y-y^2)*(z-z^2)
  void evaluate (const DomainType & x , RangeType & ret) 
  {
    ret = 1.0;
    for(int i=0; i<DomainType::dimension; i++)
      ret *= ( x[i] - SQR(x[i]) );
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) 
  {
    evaluate ( x , ret );
  }
};
 
// diffusion coefficient for this problem the id 
class Tensor : public Function < FunctionSpace < double , double, dimworld , 1 >, Tensor >
{
  typedef FunctionSpace< double , double, dimworld , 1 >  FuncSpace;
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
  typedef typename GridPartType :: IntersectionIteratorType
    IntersectionIteratorType;

  const FunctionSpaceType & space = rhs.getFunctionSpace();
  const GridPartType & gridPart = space.gridPart();
    
  typedef typename DiscreteFunctionType::DofIteratorType DofIterator;
  DofIterator dit = rhs.dbegin();
      
  IntersectionIteratorType endit = gridPart.iend(en);
  for(IntersectionIteratorType it = gridPart.ibegin(en); 
      it != endit; ++it)
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
//   // we dont not use all levels of grid for calculation, only maxlevel
// #if SGRID
//    // this leads to the same number of points for SGrid and AlbertGrid
//    int n[dimworld];
//    double h[dimworld];
//    for(int i=0; i<dimworld; i++)  { n[i] = 2; h[i] = 1.0; }

//    GridType grid ((int *) &n, (double *) &h );
// #else

   GridPtr<GridType> gridptr(filename); 
//   GridType grid ( filename );
//#endif

   gridptr->globalRefine (maxlevel);

   GridPartType part ( *gridptr );

   FuncSpaceType linFuncSpace ( part );
   std::cout << "\nSolving for " << linFuncSpace.size() << " number of unkowns. \n\n";
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
    
   bool verbose = true; 
   double dummy = 12345.67890;
   InverseOperatorType cg ( laplace, dummy , 1E-6 , 20000 , verbose );
     
   // solve linear system with cg 
   cg(rhs,solution);

   // calculation L2 error 
   ExactSolution u ( linFuncSpace ); 
   L2Error < DiscreteFunctionType > l2err;

   // pol ord for calculation the error chould by higher than 
   // pol for evaluation the basefunctions 
   double error = l2err.norm<polOrd + 2> (u ,solution, 0.0);
   std::cout << "\nL2 Error : " << error << "\n\n";

   
#if HAVE_GRAPE
   // if grape was found then display solution 
   if(turn > 0)
   {
     GrapeDataDisplay < GridType > grape(*gridptr); 
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

// #if AGRID
//   char tmp[16]; sprintf(tmp,"%d",dimworld);
//   std::string macroGridName (tmp); 
//   macroGridName += "dgrid.al";
// #else 
// #endif

//#if defined ALBERTAGRID || ALUGRID_SIMPLEX   
//  std::string macroGridName ("square_simplex.dgf");
//#else 
//  std::string macroGridName ("square_cube.dgf");
//#endif
  std::string macroGridName ("square.dgf");

  std::cout << "loading dgf " << macroGridName << "\n";
  
  ml -= refStepsForHalf;
  if(ml < 0) ml = 0;
  for(int i=0; i<2; i++)
  {
    error[i] = algorithm ( macroGridName.c_str() ,  ml , i);
    ml += refStepsForHalf ;
  }
  double eoc = log( error[0]/error[1]) / M_LN2; 
  std::cout << "EOC = " << eoc << " \n";
  return 0;
}

