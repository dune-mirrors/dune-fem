#include <config.h>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

// #include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
static const int dimw = Dune::GridSelector::dimworld;

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/misc/l2norm.hh>

#include <dune/fem/misc/double.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/solver/timeprovider.hh>

using namespace Dune;
using namespace Fem;

struct OutputParameters1
: public LocalParameter<DataOutputParameters,OutputParameters1> {
  virtual ~OutputParameters1() {}
  //! path where the data is stored (path are always relative to fem.prefix)
  virtual std::string path() const {
    return "eoc";
  }
  //! format of output (fem.io.outputformat)
  virtual int outputformat() const {
    return 1;
  }
  //! use online grape display (fem.io.grapedisplay)
  virtual bool grapedisplay() const {
    return WANT_GRAPE;
  }
  //! save data every savestep interval (fem.io.savestep)
  virtual double savestep() const {
    return -1;
  }
  //! save data every savecount calls to write method
  virtual int savecount() const {
    return 1;
  }
  //! number for first data file (no parameter available)
  virtual int startcounter() const {
    return 0;
  }
  virtual bool willWrite(bool test) const {
    return test;
  }
};

struct OutputParameters2 :
  public LocalParameter<DataOutputParameters,OutputParameters2> {
  virtual ~OutputParameters2() {}
  //! path where the data is stored (path are always relative to fem.prefix)
  virtual std::string path() const {
    return "time";
  }
  //! format of output (fem.io.outputformat)
  virtual int outputformat() const {
    return 2;
  }
  //! save data every savestep interval (fem.io.savestep)
  virtual double savestep() const {
    return 0.1;
  }
  //! save data every savecount calls to write method
  virtual int savecount() const {
    return 0;
  }
};

// polynom approximation order of quadratures,
// at least polynom order of basis functions
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
typedef GridSelector::GridType GridType;

//! the index set we are using
//typedef HierarchicGridPart<GridType> GridPartType;
//typedef DGAdaptiveLeafGridPart<GridType> GridPartType;
typedef AdaptiveLeafGridPart<GridType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
//typedef MatrixFunctionSpace < double , double, dimw , 3,5 > FuncSpace;
typedef FunctionSpace < GridType :: ctype, double , dimw , 4 > FuncSpace;

//! define the function space our unkown belong to
//! see dune/fem/lagrangebase.hh
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType,
	polOrd,CachingStorage> DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
//! dune/fem/discfuncarray.hh
typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

// the exact solution to the problem for EOC calculation
struct ExactSolution
: public Fem::Function< FuncSpace, ExactSolution >
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;

  ExactSolution ( double time = 0 )
  : time_( time )
  {}

  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 2.*cos(2.*M_PI*time_);
    for(int i=0; i<DomainType::dimension; i++)
      ret[0] *= x[i]*(1.0 -x[i]);
    for(int i=0; i<DomainType::dimension; i++)
      ret[3] *= sin(2.*M_PI*x[i]*(1.0 -x[i]));
    ret[1] = -x[1];
    ret[2] = x[0];
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }

private:
  double time_;
};

// ********************************************************************
double algorithm (GridType& grid, DiscreteFunctionType& solution,double time=0)
{
   // create exact solution for error evaluation
   ExactSolution f( time );

   // L2 error class
   L2Norm< GridPartType > l2norm( solution.space().gridPart() );

   //! perform l2-projection
   interpolate(f, solution);

   // calculation L2 error
   // pol ord for calculation the error chould by higher than
   // pol for evaluation the basefunctions
   double error = l2norm.distance( f ,solution );
   std::cout << "\nL2 Error: " << error << std::endl;
   return sqrt(error*error);
}

template <class ConsType>
struct AddLsgErr {
  typedef typename ConsType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType :: GridType ::
          template Codim<0> :: Entity EntityType;
  typedef typename EntityType :: Geometry GeometryImp;
  typedef typename ConsType :: DiscreteFunctionSpaceType :: FunctionSpaceType
                   ConsFunctionSpaceType;
  typedef typename ConsFunctionSpaceType :: DomainType ConsDomainType;
  typedef typename ConsFunctionSpaceType :: RangeType ConsRangeType;

  typedef FunctionSpace<double,double,
      ConsDomainType::dimension,
      3*ConsRangeType::dimension>
      FunctionSpaceType;
  typedef typename FunctionSpaceType :: DomainType DomainType;
  typedef typename FunctionSpaceType :: RangeType RangeType;
  typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
  /*
  virtual std::string name(int i) const {
    const int d = ConsRangeType::dimension;
    std::string ret;
    switch (i/d) {
    case 0: ret="sol-"; break;
    case 1: ret="lsg-"; break;
    case 2: ret="err-"; break;
    }
    ret+=space_.name(i%d);
    return ret;
  }
  virtual CompType comptype(int i) const {
    const int d = ConsRangeType::dimension;
    return space_.comptype(i%d);
  }
  virtual unsigned int vecsize(int i) const {
    const int d = ConsRangeType::dimension;
    return space_.vecsize(i%d);
  }
  */
  AddLsgErr(const ConsType& Uh,double time) :
    space_(Uh.space()),
    lUh_(Uh),
    time_(time),
    entity_(0),
    initialized_(false) {
  }
  ~AddLsgErr() {
  }
  template< class PointType >
  void evaluate ( const PointType &x, RangeType& ret) const {
    assert(initialized_);
    ConsRangeType u;
    lUh_.evaluate(x,u);
    ConsRangeType u0;
    DomainType global = entity().geometry().global( coordinate( x ) );
    ExactSolution f( time_ );
    f.evaluate(global,u0);
    const int d = ConsRangeType::dimension;
    for (int i=0;i<d;i++) {
      ret[0*d+i] = u[i];
      ret[1*d+i] = u0[i];
      ret[2*d+i] = std::abs(u[i]-u0[i]);
    }
  }
  template< class PointType >
  void jacobian ( const PointType &x, JacobianRangeType& ret) const {
    abort();
  }
  //! init local function
  void init(const EntityType& en)
  {
    entity_ = &en;
    lUh_.init(en);
    initialized_ = true;
  }
private:
  const EntityType &entity() const
  {
    assert( entity_ );
    return *entity_;
  }

  typedef typename ConsType::LocalFunctionType ConsLocalFunctionType;
  const DiscreteFunctionSpaceType& space_;
  ConsLocalFunctionType lUh_;
  double time_;
  const EntityType* entity_;
  bool initialized_;
};

/*
struct Model { // : public SpaceDescriptor {
  virtual std::string name(int i) const {
    switch (i) {
    case 0: return "rho";
    case 1:
    case 2: return "momentum";
    case 3: return "energy";
    }
  }
  virtual CompType comptype(int i) const {
    switch (i) {
    case 0: return scalar;
    case 1:
    case 2: return vector;
    case 3: return scalar;
    }
  }
  virtual unsigned int vecsize(int i) const {
    switch (i) {
    case 0: return 0;
    case 1:
    case 2: return 2;
    case 3: return 0;
    }
  }
};

*/
//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main (int argc, char **argv)
{
  MPIManager :: initialize( argc, argv );
  Parameter::append(argc,argv);
  try
  {
    int ml = (argc > 2)? atoi( argv[1] ) : 2;
    Parameter::append("parameter");
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
    // Model model;
    // solution.space().setDescription(model);

    typedef AddLsgErr<DiscreteFunctionType> AddLsgErrType;
    AddLsgErrType evalAddLsgErr(solution,0);

    typedef LocalFunctionAdapter<AddLsgErrType> AddLsgErrFunction;
    AddLsgErrFunction addLsgErr("U",evalAddLsgErr,solution.space().gridPart());

    {
      auto output = dataOutput( grid, addLsgErr ); //, OutputParameters1() );
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
      // tp.setEndTime(1);
      auto output = dataOutput( grid, tp, addLsgErr, OutputParameters2() );
      for( tp.init(0.01) ; tp.time()<=0.04 ; tp.next(0.01) )
      {
        algorithm ( grid , solution,tp.time() );
        output.write(tp);
      }
      output.write();
    }

    return 0;
  }
  catch( const Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

