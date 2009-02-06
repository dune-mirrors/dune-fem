#include <iostream>
#include <config.h>
#include <string>
#include <sstream>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
static const int dimw = dimworld;

#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh> 
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh> 
#include <dune/fem/space/common/adaptmanager.hh>

#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/misc/l2error.hh>

#include <dune/fem/misc/double.hh>

#include <dune/fem/io/file/dataoutput.hh>

using namespace Dune;
struct OutputParameters1 : public DataOutputParameters {
  virtual ~OutputParameters1() {}
  //! path where the data is stored (path are always relative to fem.prefix)
  virtual std::string path() const {
    return "data";
  }
  //! base of file name for data file (fem.io.datafileprefix)
  virtual std::string prefix() const {
    return "testEOC";
  }
  //! format of output (fem.io.outputformat)
  virtual int outputformat() const {
    return 2;
  }
  //! use online grape display (fem.io.grapedisplay)
  virtual bool grapedisplay() const {
    return WANT_GRAPE;
  }
  //! save data every savestep interval (fem.io.savestep)
  virtual double savestep() const {
    return -1;
  }
  //! save data every savecount calls to write method (fem.io.savecount)
  virtual int savecount() const {
    return 1;
  }
  //! number for first data file (no parameter available)
  virtual int startcounter() const {
    return 0;
  }
};
struct OutputParameters2 : public DataOutputParameters {
  virtual ~OutputParameters2() {}
  //! path where the data is stored (path are always relative to fem.prefix)
  virtual std::string path() const {
    return "data";
  }
  //! base of file name for data file (fem.io.datafileprefix)
  virtual std::string prefix() const {
    return "testTIME";
  }
  //! format of output (fem.io.outputformat)
  virtual int outputformat() const {
    return 2;
  }
  //! use online grape display (fem.io.grapedisplay)
  virtual bool grapedisplay() const {
    return WANT_GRAPE;
  }
  //! save data every savestep interval (fem.io.savestep)
  virtual double savestep() const {
    return 0.1;
  }
  //! save data every savecount calls to write method (fem.io.savecount)
  virtual int savecount() const {
    return -1;
  }
  //! number for first data file (no parameter available)
  virtual int startcounter() const {
    return 0;
  }
};

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
//typedef HierarchicGridPart<GridType> GridPartType;
//typedef DGAdaptiveLeafGridPart<GridType> GridPartType;
typedef AdaptiveLeafGridPart<GridType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
//typedef MatrixFunctionSpace < double , double, dimw , 3,5 > FuncSpace;
typedef FunctionSpace < GridType :: ctype, double , dimw , 2 > FuncSpace;

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
  double time_;
public:
  ExactSolution (const FuncSpace &f,double time=0) : 
    Function < FuncSpace , ExactSolution > ( f ),
    time_(time) {}
 
  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 2.*cos(2.*M_PI*time_); // maximum of function is 2
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
double algorithm (GridType& grid, DiscreteFunctionType& solution,double time=0)
{
   // create exact solution for error evaluation 
   ExactSolution f ( solution.space(),time ); 

   // L2 error class 
   L2Error < DiscreteFunctionType > l2err;
       
   //! perform l2-projection
   L2ProjectionImpl::project(f, solution);

   // calculation L2 error 
   // pol ord for calculation the error chould by higher than 
   // pol for evaluation the basefunctions 
   typedef DiscreteFunctionSpaceType :: RangeType RangeType; 
   RangeType error = l2err.norm(f ,solution, 0.0);
   for(int i=0; i<RangeType::dimension; ++i)
     std::cout << "\nL2 Error["<<i<<"] : " << error[i] << "\n\n";

   return sqrt(error*error);
}


//**************************************************
//
//  main programm, run algorithm twice to calc EOC 
//
//**************************************************
int main (int argc, char **argv)
{
  MPIManager :: initialize( argc, argv );
  try
  {
  
  if(argc != 2)
  {
    fprintf(stderr,"usage: %s <maxlevel> \n",argv[0]);
    exit(1);
  }
  int ml = atoi( argv[1] );
  std::vector< double> error(ml);
  std::stringstream tmp; 
  tmp << dimw;
  std::string macroGridName (tmp.str()); 
  macroGridName += "dgrid.dgf";

  GridPtr<GridType> gridptr(macroGridName);
  GridType& grid=*gridptr;
  const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();

  GridPartType part ( grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  solution.clear();
  
  typedef Tuple<DiscreteFunctionType*> OutputType;
  OutputType out(&solution);
  {
    DataOutput<GridType,OutputType> output(grid,out,OutputParameters1());
    for(int i=0; i<ml; i+=step)
    {
      GlobalRefine::apply(grid,step);
      error[i] = algorithm ( grid , solution );
      output.write();
      if (i>0) 
      {
        double eoc = log( error[i-step]/error[i]) / M_LN2; 
        std::cout << "EOC = " << eoc << " \n";
      }
    }
  }
  {
    GridTimeProvider<GridType> tp(0,grid);
    tp.setEndTime(1);
    DataOutput<GridType,OutputType> output(grid,out,tp,OutputParameters2());
    for( tp.init(0.01) ; tp ; tp.next(0.01) )
    {
      algorithm ( grid , solution,tp.time() );
      output.write(tp);
    }
    output.writeData();
  }
  return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

