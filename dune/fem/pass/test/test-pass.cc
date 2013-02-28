#include <config.h>

#include <iostream>
#include <cstdlib>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolation.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/test/exactsolution.hh>
#include <dune/fem/test/testgrid.hh>

#include <dune/fem/pass/applylocaloperator.hh>
#include <dune/fem/pass/common/wrapper.hh>

// Internal forward declaration
// ----------------------------

template< class, int >
class DiscreteModel;



// DiscreteModelTraits
// -------------------

template< class DiscreteFunction, int id >
struct DiscreteModelTraits
{
  typedef DiscreteModel< DiscreteFunction, id > DiscreteModelType;

  typedef DiscreteFunction DestinationType;

  typedef typename DestinationType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
};



// DiscreteModel
// -------------

template< class DiscreteFunction, int id >
struct DiscreteModel
: public Dune::Fem::ApplyLocalOperatorDiscreteModel< DiscreteModelTraits< DiscreteFunction, id >, id >
{
  typedef DiscreteModelTraits< DiscreteFunction, id > Traits;

private:
  typedef Dune::Fem::ApplyLocalOperatorDiscreteModel< Traits, id > BaseType;
  
  Dune::integral_constant< int, id > u;

public:
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalCoordinateType LocalCoordinateType;

  DiscreteModel ( int order ) : order_( order ), time_( 0. ) {}

  void setTime ( double time ) { time_ = time; }

  double time () const { return time_; }

  int order ( const EntityType &entity ) const { return order_; }

  template< class RangeTuple >
  void evaluate ( const EntityType &entity,
                  const LocalCoordinateType &x,
                  const RangeTuple &tuple,
                  RangeType &value ) const
  {
    Dune::Fem::Wrapper< RangeTuple, typename BaseType::Selector > wrapper( tuple );
    value = wrapper[ u ];
  }

  template< class JacobianTuple >
  void jacobian ( const EntityType &entity,
                  const LocalCoordinateType &x,
                  const JacobianTuple &tuple,
                  JacobianRangeType &jacobian ) const
  {
    Dune::Fem::Wrapper< JacobianTuple, typename BaseType::Selector > wrapper( tuple );
    jacobian = wrapper[ u ];
  }

private:
  int order_;
  double time_;
};



// LocalRestrictionOperator
// ------------------------

template< class DiscreteFunctionSpace >
struct LocalRestrictionOperator
{
  typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

private:
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType > FiniteVolumeSpaceType;
  typedef Dune::Fem::TemporaryLocalFunction< FiniteVolumeSpaceType > TemporaryLocalFunctionType;

public:
  LocalRestrictionOperator ( const DiscreteFunctionSpaceType &space )
  : space_( space ),
    fvSpace_( space.gridPart() ),
    temporary_( fvSpace_ )
  {}

  template< class LocalFunction, class LocalDofVector >
  void operator() ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
  {
    const typename DiscreteFunctionSpaceType::EntityType &entity = localFunction.entity();
    temporary_.init( entity );
    fvSpace_.interpolate( localFunction, temporary_ );
    space_.interpolate( temporary_, dofs );
  }

private:
  const DiscreteFunctionSpaceType &space_;
  FiniteVolumeSpaceType fvSpace_;
  mutable TemporaryLocalFunctionType temporary_;
};



// Operator
// --------

template< class DiscreteFunction >
class Operator
: public Dune::Fem::SpaceOperatorInterface< DiscreteFunction >
{
  typedef Dune::Fem::SpaceOperatorInterface< DiscreteFunction > BaseType;

public:
  typedef typename BaseType::DestinationType DiscreteFunctionType;
  typedef typename BaseType::SpaceType DiscreteFunctionSpaceType;

private:
  enum PassId{ Start, Interpolation };

  typedef DiscreteModel< DiscreteFunctionType, Start > DiscreteModelType;
  typedef LocalRestrictionOperator< DiscreteFunctionSpaceType > LocalOperatorType;

  typedef Dune::Fem::StartPass< DiscreteFunctionType, Start > StartPassType;
  typedef Dune::Fem::ApplyLocalOperatorPass< DiscreteModelType, LocalOperatorType, StartPassType, Interpolation > InterpolationPassType;

public:
  Operator ( const DiscreteFunctionSpaceType &space )
  : space_( space ),
    discreteModel_( space.order() ),
    localOperator_( space ),
    interpolationPass_( discreteModel_, localOperator_, startPass_, space_ )
  {}

  void operator () ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest ) const
  {
    interpolationPass_( arg, dest );
  }

  void setTime( double time ) { interpolationPass_.setTime( time ); }

  double timeStepEstimate () const { return interpolationPass_.timeStepEstimate(); }

  const DiscreteFunctionSpaceType &space () const { return space_; }

private:
  const DiscreteFunctionSpaceType &space_;
  DiscreteModelType discreteModel_;
  LocalOperatorType localOperator_;
  StartPassType startPass_;
  InterpolationPassType interpolationPass_;
};



// Typedefs used in the following
// ------------------------------

typedef Dune::GridSelector::GridType GridType;

const int dimRange = DIMRANGE;
typedef Dune::Fem::FunctionSpace< GridType::ctype, double, GridType::dimensionworld, dimRange > FunctionSpaceType;
typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;

const int order = POLORDER;
typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, order > DGDiscreteFunctionSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< DGDiscreteFunctionSpaceType > DGDiscreteFunctionType;

typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType > FVDiscreteFunctionSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< FVDiscreteFunctionSpaceType > FVDiscreteFunctionType;



int main ( int argc, char **argv )
try
{
  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameter file
  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  // get grid
  GridType &grid = Dune::Fem::TestGrid::grid();

  // initial refinement
  int refCount = Dune::Fem::Parameter::getValue< int >( "startlevel", 0 );
  const int refineStepsForHalf = Dune::Fem::TestGrid::refineStepsForHalf();
  Dune::Fem::GlobalRefine::apply( grid, refCount*refineStepsForHalf );

  // create grid part
  GridPartType gridPart( grid );
  
  // discrete function spaces
  DGDiscreteFunctionSpaceType dgSpace( gridPart );
  FVDiscreteFunctionSpaceType fvSpace( gridPart );
 
  // create grid function
  typedef Dune::Fem::ExactSolution< FunctionSpaceType > FunctionType;
  FunctionType function;
  typedef Dune::Fem::GridFunctionAdapter< FunctionType, GridPartType > GridFunctionType;
  GridFunctionType gridFunction( "Exact solution", function, gridPart, 2*order );

  // create discrete functions
  FVDiscreteFunctionType u( "FV function u", fvSpace );
  DGDiscreteFunctionType v( "DG function v", dgSpace );
  DGDiscreteFunctionType w( "DG function w", dgSpace );

  // compute projections
  Dune::Fem::Interpolation< FVDiscreteFunctionType >::apply( gridFunction, u );
  Dune::Fem::Interpolation< DGDiscreteFunctionType >::apply( gridFunction, v );

  // apply operator
  Operator< DGDiscreteFunctionType > op( dgSpace );
  op( v, w );

  // write output
  typedef Dune::tuple< const DGDiscreteFunctionType *, const DGDiscreteFunctionType * > IOTupleType;
  IOTupleType ioTuple( &v, &w );
  typedef Dune::Fem::DataOutput< GridType, IOTupleType > DataOutputType;
  DataOutputType dataOutput( grid, ioTuple );
  dataOutput.write();

  // compute error
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  double error = l2norm.distance( u, w ); 
  std::cout << "Error ||u - w||_L2 = " << error << std::endl;

  return 0;
}
catch( Dune::Exception &e )
{
  std::cerr << "Dune exception thrown: " << e.what() << std::endl;
  return 1;
}
catch ( std::exception &e )
{
  std::cerr << "Exception thrown: " << e.what() << std::endl;
  return 1;
}
catch (...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}

