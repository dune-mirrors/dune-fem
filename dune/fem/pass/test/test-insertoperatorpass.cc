#include <config.h>

#include <cstdlib>

#include <array>
#include <iostream>
#include <memory>
#include <utility>

#include <dune/common/fvector.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#include <dune/fem/pass/common/pass.hh>
#include <dune/fem/pass/insertoperator.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/test/exactsolution.hh>

// OperatorWrapper
// ---------------

template< class Operator >
class OperatorWrapper final
  : public Dune::Fem::Operator< typename Operator::DomainFunctionType, typename Operator::RangeFunctionType >
{
  typedef Dune::Fem::Operator< typename Operator::DomainFunctionType, typename Operator::RangeFunctionType > BaseType;

public:
  typedef typename BaseType::DomainFunctionType DomainFunctionType;
  typedef typename BaseType::RangeFunctionType RangeFunctionType;

private:
  enum { start, insertOperator };

  typedef Dune::Fem::StartPass< DomainFunctionType, start > StartPassType;
  typedef Dune::Fem::InsertOperatorPass< Operator, StartPassType, insertOperator > InsertOperatorPassType;

public:
  OperatorWrapper ( const typename RangeFunctionType::DiscreteFunctionSpaceType &space,
                    const Operator &op )
    : insertOperatorPass_( space, op, startPass_ )
  {}

  void operator() ( const DomainFunctionType &u, RangeFunctionType &v ) const override
  {
    insertOperatorPass_( u, v );
  }

private:
  StartPassType startPass_;
  InsertOperatorPassType insertOperatorPass_;
};



// project
// -------

template< class DomainFunction, class RangeFunction >
void project ( const DomainFunction &vh, RangeFunction &wh )
{
  typedef Dune::Fem::VtxProjection< DomainFunction, RangeFunction > VertexProjectionType;
  VertexProjectionType projection;

  typedef OperatorWrapper< VertexProjectionType > WrappedOperatorType;
  WrappedOperatorType wrapped( wh.space(), projection );
  wrapped( vh, wh );
}



// createGrid
// ----------

template< class ctype, int dim >
std::unique_ptr< Dune::YaspGrid< dim, Dune::EquidistantOffsetCoordinates< ctype, dim > > >
createGrid( int count = 1 )
{
  // select unit cube as domain
  typedef Dune::FieldVector< ctype, dim > GlobalCoordinate;
  GlobalCoordinate a( 0 ), b( 1 );

  // macro level has one element
  std::array< int, dim > cells;
  cells.fill( 1 );

  // create and refine grid
  typedef Dune::YaspGrid< dim, Dune::EquidistantOffsetCoordinates< ctype, dim > > GridType;
  auto grid = std::make_unique< GridType >( a, b, cells );
  grid->loadBalance();
  grid->globalRefine( count );
  grid->loadBalance();

  return grid;
}



int main ( int argc, char **argv )
{
  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameters
  int count = 1;
  if( argc > 1 )
    count = std::atoi( argv[ 1 ] );

  // create grid
  auto grid = createGrid< double, 2 >( count );
  typedef decltype( grid )::element_type GridType;

  // choose grid part
  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( *grid );

  // create analytical function
  typedef Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 3 > FunctionSpaceType;
  typedef Dune::Fem::ExactSolution< FunctionSpaceType > FunctionType;
  FunctionType function;

  // adapter to grid function
  typedef Dune::Fem::GridFunctionAdapter< FunctionType, GridPartType > GridFunctionType;
  GridFunctionType uh( "uh", function, gridPart, 2 );

  // piecewise constant domain function
  typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType > DomainFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DomainFunctionSpaceType > DomainFunctionType;
  DomainFunctionSpaceType domainSpace( gridPart );
  DomainFunctionType vh( "vh", domainSpace );

  // continuous range function
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 > RangeFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< RangeFunctionSpaceType > RangeFunctionType;
  RangeFunctionSpaceType rangeSpace( gridPart );
  RangeFunctionType wh( "wh", rangeSpace );

  // initialize domain function
  Dune::Fem::interpolate( uh, vh );
  project( vh, wh );

  return 0;
}
