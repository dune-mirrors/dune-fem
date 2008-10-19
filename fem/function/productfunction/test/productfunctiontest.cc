#include <iostream>
#include <config.h>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
static const int dimw = dimworld;
static const int dimp = dimworld;
  
//dimension of second space/ grid
static const int dim2 = 3;

//polynom order of second space
// static const int polOrd2 = 4;
//#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#include <dune/grid/io/file/dgfparser/dgfs.hh>

#include <dune/fem/function/productfunction.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh> 

#include <dune/grid/common/referenceelements.hh>
 
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif  
using namespace Dune;
 
// polynom approximation order of quadratures, 
// at least polynom order of basis functions  
const int polOrd = POLORDER;

//***********************************************************************
/*! L2 Projection of a function f on XxY: 
  The Projection should converge to the given function f.
  with the finite element method using lagrangian elements of polynom order +1.
*/
//***********************************************************************

// static const int dimworld = GRIDDIM;
// static const int dimworld = dimworld;
typedef SGrid  < dim2,dim2 > Grid2Type;
// static const int refStepsForHalf = 1;
//! the index set we are using 
typedef HierarchicGridPart<GridType> GridPartType;
typedef HierarchicGridPart<Grid2Type> GridPart2Type;
 
//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
typedef FunctionSpace < double , double, dimp , 1 > FuncSpace;
typedef FunctionSpace < double , double, dim2 , 1 > FuncSpace2;
//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh

//! Get the Dofmanager type
typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType; 

template <int polOrd2> 
struct TestProdDF {
  typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
				     polOrd,CachingStorage> DiscreteFunction1SpaceType;
  /* 
     typedef DiscontinuousGalerkinSpace<FuncSpace2, GridPart2Type, 
     polOrd2 ,CachingStorage> DiscreteFunction2SpaceType; 
  */
  typedef LegendreDiscontinuousGalerkinSpace<FuncSpace2, GridPart2Type, 
					     polOrd2 ,CachingStorage> DiscreteFunction2SpaceType; 
  //! define the type of discrete function we are using , see
  //! dune/fem/discfuncarray.hh
  typedef ProductDiscreteFunction < DiscreteFunction1SpaceType, DiscreteFunction2SpaceType > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunction1Type  DiscreteFunction1Type;
  //! the exact solution to the problem for EOC calculation 
  class ExactSolution : public Function < FuncSpace , ExactSolution > 
  {   
    typedef typename FuncSpace::RangeType RangeType; 
    typedef typename FuncSpace::RangeFieldType RangeFieldType;
    typedef typename FuncSpace::DomainType DomainType;     
    typedef typename FuncSpace2::DomainType Domain2Type;
  public: 
    ExactSolution (FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {} 
  
    //! f(x,y) = 4x*(1-x)*y*(1-y)  
    void evaluate (const DomainType & x , RangeType & ret)  const 
    {
      ret =1.; 
   
      for(int i=0; i<DomainType::dimension; i++) 
      ret *= x[i]*(1.0 -x[i])*4.0;
      //ret *= sin( x[i]);
      //ret += x[i];
    
    } 
  
    void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  
    void setPara (const DomainType & x ,const Domain2Type & y , RangeType & ret)  const
    {
      evaluate ( x , ret ); 
      for(int i=0; i<Domain2Type::dimension; i++)
	{ 
	  ret *=1.5*cos(1.33 *M_PI*y[i]);   //! Note: BaseFunctions must be ortho-normal!!!! 
	  //ret *=1.0- 0.2* y[i];
	  //ret *= pow(y[i],4);
	}   
      //std::cerr << x << " " << y << " " << ret << std::endl;
    }
  
  };   
    

  
  // ********************************************************************
  template <class DiscreteFunctionType>
  class L2Projection
  {
    typedef typename DiscreteFunction1Type::FunctionSpaceType DiscreteFunction1SpaceType;
    typedef typename DiscreteFunctionType::DiscreteFunction1Type  DiscreteFunction1Type;

  public:
    template <class FunctionType>
    static void project (const FunctionType &f, DiscreteFunctionType &discFunc, int polOrd) 
    {
      typedef typename DiscreteFunction1SpaceType::Traits::GridType GridType;
      typedef typename DiscreteFunction1SpaceType::Traits::IteratorType Iterator;
    
      typedef typename DiscreteFunction2SpaceType::Traits::GridType Grid2Type;
      typedef typename DiscreteFunction2SpaceType::Traits::IteratorType Iterator2;
    
      const DiscreteFunction1SpaceType& space =  discFunc.space();
      const DiscreteFunction2SpaceType& space2 =  discFunc.space2();
    
      discFunc.clear(); 
    
      //! Note: BaseFunctions must be ortho-normal!!!! 
      typedef typename DiscreteFunction1SpaceType::BaseFunctionSetType BaseFunctionSetType ; 
      typedef typename DiscreteFunction2SpaceType::BaseFunctionSetType BaseFunctionSet2Type;

      typedef typename DiscreteFunction1Type::LocalFunctionType LocalFuncType;

      typename DiscreteFunction1SpaceType::RangeType ret (0.0);
      typename DiscreteFunction1SpaceType::RangeType phi (0.0);
      typename DiscreteFunction2SpaceType::RangeType psi (0.0);
    
      int quadOrd2 = 2*space2.order()+1;
    
      Iterator endit = space.end();
      Iterator2 endit2 = space2.end(); 
    
      for(Iterator it = space.begin(); it != endit ; ++it) 
	{
	  ElementQuadrature<GridPartType,0> quad(*it, polOrd);
	  // CachingQuadrature<GridPartType,0> quad(*it, polOrd);
	  const int quadNop  = quad.nop();
	  for(Iterator2 it2 = space2.begin(); it2 != endit2 ; ++it2) 
	    {
	      // Get quadrature rule
	      ElementQuadrature<GridPart2Type,0> quad2(*it2, quadOrd2);
	      // CachingQuadrature<GridPart2Type,0> quad2(*it2, quadOrd2);
      
	      const BaseFunctionSet2Type bSet2 = space2.baseFunctionSet(*it2);
	      const int numOfDofs = bSet2.numBaseFunctions();
	      const int quadNop2 = quad2.nop();
	      const typename Grid2Type::template Codim<0>::Entity::Geometry& itGeom2 = 
		(*it2).geometry();
	
	      for(int j=0; j< numOfDofs; j++) 
		{
		  DiscreteFunction1Type ldf = discFunc.localFunction(*it2,j);  
        
		  // Get quadrature rule
		  // const int quadNop  = quad.nop();
		  LocalFuncType lf = ldf.localFunction(*it);  
        
		  const BaseFunctionSetType& baseset = lf.baseFunctionSet();
	
		  const typename GridType::template Codim<0>::Entity::Geometry& 
		    itGeom = (*it).geometry();
	
		  const int numDofs = lf.numDofs();
	
		  for(int qP2 = 0; qP2 < quadNop2 ; ++qP2) 
		    {
		      bSet2.evaluate(j,quad2[qP2],psi);
	    
		      for(int qP = 0; qP < quadNop ; ++qP) 
			{
			  f.setPara(itGeom.global(quad.point(qP))
				    ,itGeom2.global(quad2.point(qP2)), ret); 
		
			  for(int i=0; i<numDofs; ++i) 
			    {
			      baseset.evaluate(i,quad[qP],phi);
			      lf[i] +=quad2.weight(qP2) * quad.weight(qP) * (ret * phi) * psi ;
			      /*std::cerr << "Project: " 
				<< i << " " << j << " " 
				<< qP << " " << qP2 << " " 
				<< phi << " " << psi << " " << ret << "  " 
				<< lf[i] << std::endl;
			      */		 
			    }
			}
		    } 
		}
	    }
	} 
    }
  

    template <class FunctionType>
    static void project (const FunctionType &f, DiscreteFunctionType &discFunc) 
    {
      const DiscreteFunction1SpaceType& space =  discFunc.space();
      int polOrd = 2 * space.order() ;
      project(f,discFunc,polOrd);
    }
  };


  // calculates || u-u_h ||_L2
  template <class DiscreteFunctionType>
  class L2Error
  {
    typedef typename DiscreteFunctionType::DiscreteFunction1Type DiscreteFunction1Type;
    typedef typename DiscreteFunction1Type::FunctionSpaceType DiscreteFunction1SpaceType;
    typedef typename DiscreteFunction1SpaceType :: RangeType RangeType;

  public:
    template <class FunctionType>
    RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
		    double time, int polOrd) const
    {
      const DiscreteFunction1SpaceType & space = discFunc.space();
      const DiscreteFunction2SpaceType& space2 =  discFunc.space2();

      typedef typename DiscreteFunction1SpaceType::GridType GridType;
      typedef typename DiscreteFunction1SpaceType::IteratorType IteratorType;
      typedef typename DiscreteFunction1Type::LocalFunctionType LocalFuncType;
      typedef typename DiscreteFunction2SpaceType::IteratorType Iterator2Type;
      typedef typename DiscreteFunction2SpaceType::BaseFunctionSetType BaseFunctionSet2Type;
    
      RangeType ret (0.0); 
      RangeType phi (0.0);
      RangeType error(0.0);

      enum { dimRange = DiscreteFunction1SpaceType :: DimRange };
      int quadOrd2 = 2*space2.order()+2;

      IteratorType endit = space.end();
      Iterator2Type endit2 = space2.end();  
      DiscreteFunction1Type ldf("tmp",space);   
    
      for(IteratorType it = space.begin(); it != endit ; ++it)
	{ 
	  // Get quadrature rule
	  ElementQuadrature<GridPartType,0> quad(*it, polOrd);
	  // CachingQuadrature<GridPartType,0> quad(*it,polOrd);
	  const int quadNop = quad.nop();
	  for(Iterator2Type it2 = space2.begin(); it2 != endit2 ; ++it2) 
	    {
	      // Get quadrature rule
	      ElementQuadrature<GridPart2Type,0> quad2(*it2,quadOrd2);
	      // CachingQuadrature<GridPart2Type,0> quad2(*it2,quadOrd2);
	      const int quadNop2 = quad2.nop();  
	  
	      for(int qP2 = 0; qP2 < quadNop2 ; ++qP2) 
		{
		  //discFunc.localFunction(*it2,quad2.point(qP2),ldf);
		  discFunc.localFunction(*it2,quad2,qP2,ldf);
		  LocalFuncType lf = ldf.localFunction(*it);
	
		  /*       #if HAVE_GRAPE
		  // if Grape was found, then display last solution 
		  std::cout << " np = " << qP2 ;
		  std::cout << " y_loc = " << quad2.point(qP2) ;
		  std::cout << " y = " << (*it2).geometry().global(quad2.point(qP2)) << std::endl;
	
		  GrapeDataDisplay < GridType > grape(space.grid()); 
		  grape.dataDisplay( ldf );
     
		  #endif 
		  */	
 
		  for(int qP = 0; qP < quadNop; ++qP)
		    {
		      double weight = 
			quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP))*quad2.weight(qP2) * (*it2).geometry().integrationElement(quad2.point(qP2)); 

		      lf.evaluate(quad[qP],phi);
   
		      f.setPara((*it).geometry().global(quad.point(qP)),(*it2).geometry().global(quad2.point(qP2)), ret);
		      /*std::cerr << "error : "
			<< qP << " " << qP2 << " "
			<< phi[0] << " " << ret[0] << " "
			<< std::endl;
		      */		    
		      for(int i=0; i< dimRange; ++i)
			{ 
			  error[i] += weight * SQR(ret[i] - phi[i]);
			}
		    }  
		}
	    }
	}
      /*#if HAVE_GRAPE
      // if Grape was found, then display last solution 
      std::cout << " np = " << qP2 ;
      std::cout << " y_loc = " << quad2.point(qP2) ;
      std::cout << " y = " << (*it2).geometry().global(quad2.point(qP2)) << std::endl;
	
	
	
      GrapeDataDisplay < GridType > grape(space.grid()); 
      grape.dataDisplay( ldf );
     
      #endif
      */
      
    

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
      const DiscreteFunction1SpaceType & space = discFunc.space();
      int polOrd = 2 * space.order() + 2;
      //std::cout << "polord = " << polOrd << " \n\n";
      return norm(f,discFunc,time,polOrd);
    }
  };
  // ********************************************************************
  double algorithm (GridType& grid, Grid2Type& grid2 ,DiscreteFunctionType& solution  , int turn )
  {
    GridPartType part ( grid );
    DiscreteFunction1SpaceType linFuncSpace ( part ); 
    GridPart2Type part2 ( grid2 );
    DiscreteFunction2SpaceType linFuncSpace2 ( part2 );
   
   
    ExactSolution f ( linFuncSpace ); 
    L2Error < DiscreteFunctionType > l2err;
  
    // calculation L2 error 
    // pol ord for calculation the error should by higher than 
    // pol for evaluation the basefunctions 
    typedef DiscreteFunction1SpaceType :: RangeType RangeType; 
   
    //! perform l2-projection
    L2Projection<DiscreteFunctionType>::
      project(f,solution ); 
 
    RangeType error = l2err.norm(f ,solution, 0.0);

    // for(int i=0; i<RangeType::dimension; ++i)
    //   std::cout << "\nL2 Error["<<i<<"] : " << error[i] << "\n\n";
 
    return error;
  }
  //**************************************************
  //
  //  main programm, run algorithm twice to calc EOC 
  //
  //**************************************************
  double main (GridType& grid, Grid2Type& grid2) {
    try {
      GridPartType part ( grid );
      GridPart2Type part2 ( grid2 );
      DiscreteFunction1SpaceType linFuncSpace ( part ); 
      DiscreteFunction2SpaceType linFuncSpace2 ( part2 );
      double error;
      
      DiscreteFunctionType solution ( linFuncSpace, linFuncSpace2 );
      solution.clear();
      // DofManagerType& dm = DofManagerFactoryType :: getDofManager( grid );
      // dm.resize();
      error = algorithm ( grid ,grid2, solution , 0);
      return error;
    }
    catch (std::exception & e) {
      std::cout << "STL ERROR: " << e.what() << std::endl;
      return 1;
    }
    catch (Dune::Exception & e) {
      std::cout << "DUNE ERROR: " << e.what() << std::endl;
      return 1;
    }
    catch (...) {
      std::cout << "Unknown ERROR" << std::endl;
      return 1;
      
    }
    return 1;
  }
};
int main (int argc, char **argv)
{
  if(argc < 2)
  {
    std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
    return 1;
  }
  const int ml = atoi( argv[1] );

  double* error = new double[8*ml];
  char tmp[16]; sprintf(tmp,"%d",dimp);
  std::string macroGridName (tmp); 
  macroGridName += "dgrid.dgf"; 
  const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();   //step=1
  GridPtr<GridType> gridptr(macroGridName);
  GridType& grid=*gridptr;
  int *N = new int[dim2];
  sgrid_ctype * H = new sgrid_ctype[dim2];
  for(int i=0; i<dim2; i++) {
    N[i]= 1;
    H[i]= 1.;
  }  
  Grid2Type grid2(N, H);
  delete [] H;
  delete [] N;
  
  for(int i=0; i<ml; i+=step) {
    grid.globalRefine(step);
    TestProdDF<1> test1;
    error[6*i+0] = test1.main(grid,grid2);
    TestProdDF<2> test2;
    error[6*i+1] = test2.main(grid,grid2);
    TestProdDF<3> test3;
    error[6*i+2] = test3.main(grid,grid2);
    TestProdDF<4> test4;
    error[6*i+3] = test4.main(grid,grid2);
    TestProdDF<5> test5;
    error[6*i+4] = test5.main(grid,grid2);
    TestProdDF<6> test6;
    error[6*i+5] = test6.main(grid,grid2);
    TestProdDF<7> test7;
    error[6*i+6] = test7.main(grid,grid2);
    TestProdDF<8> test8;
    error[6*i+7] = test8.main(grid,grid2);
    /*
      if (i>0) {
      double eoc = log( error[i-step]/error[i]) / M_LN2;  
      //double eoc2 =  error[i-step]/error[i];
      std::cout << "EOC = " << eoc << " \n\n";
      //std::cout << "fehler alt:neu = " << eoc2 << " \n\n";
      
    }
    */
    for (int j=0;j<8;j++) 
      std::cout << i << " " << j << " " << error[6*i+j] << " # Result" 
		<< std::endl;
    std::cout << std::endl;
  }
}
