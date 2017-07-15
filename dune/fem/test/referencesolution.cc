#include <config.h>

#include <dune/fem/gridpart/levelgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/fem/test/testgrid.hh>
#include <dune/fem/test/exactsolution.hh>

template< class GF, class GP >
struct ToGridPart : public Dune::Fem::BindableFunction< GP, typename GF::RangeType >
{
  typedef typename GF::RangeType Range;
  typedef Dune::Fem::BindableFunction< GP, Range > Base;
  typedef typename GF::GridPartType OtherGridPartType;
  typedef typename Base::EntityType EntityType;

  ToGridPart ( const GF &gf )
    : lgf_(gf), gridPart_(gf.gridPart())
  {}

  template< class Point >
  void evaluate ( const Point &p, typename Base::RangeType &val ) const
  {
    auto x = Dune::Fem::coordinate(p);
    EntityType coarse = fine_;
    while ( !gridPart_.contains(coarse) )
    {
      x = coarse.geometryInFather().global(x);
      coarse = coarse.father();
    }
    lgf_.evaluate(x,val);
  }
  template< class Point >
  void evaluate ( const Point &p, typename Base::JacobianRangeType &val ) const
  {
    auto x = Dune::Fem::coordinate(p);
    EntityType coarse = fine_;
    while ( !gridPart_.contains(coarse) )
    {
      x = coarse.geometryInFather().global(x);
      coarse = coarse.father();
    }
    lgf_.jacobian(x,val);
  }
  void bind ( const EntityType &entity )
  {
    EntityType coarse = entity;
    fine_ = entity;
    while ( !gridPart_.contains(coarse) ) // missing
      coarse = coarse.father();
    lgf_.bind( coarse );
  }
  void unbind() { lgf_.unbind(); }
  unsigned int order() const { return lgf_.order(); }

private:
  Dune::Fem::ConstLocalFunction< GF > lgf_;
  const OtherGridPartType &gridPart_;
  EntityType fine_;
};

int main ( int argc, char ** argv )
try
{
  typedef Dune::GridSelector::GridType GridType;

  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

  unsigned int maxLevel = 5;
  auto& grid = Dune::Fem::TestGrid::grid();
  grid.globalRefine( Dune::DGFGridInfo< GridType >::refineStepsForHalf()*maxLevel );

  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );
  typedef Dune::Fem::FunctionSpace< GridType::ctype, GridType::ctype, GridType::dimensionworld, 1 > FunctionSpaceType;
#if HAVE_DUNE_LOCALFUNCTIONS
  typedef Dune::Fem::LagrangeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
  const int polOrder = Dune::Fem::Parameter::getValue< int >( "fem.lagrange.polinomialOrder");
  DiscreteFunctionSpaceType space( gridPart, polOrder );
#else
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2 > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType space( gridPart );
#endif
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > FunctionType;
  FunctionType reference( "reference", space );
  interpolate( gridFunctionAdapter( Dune::Fem::ExactSolution<FunctionSpaceType>(), gridPart, space.order() + 2 ), reference );

  std::vector< Dune::FieldVector<double,2> > error;
  for (unsigned int level=0;level<=maxLevel;++level)
  {
    typedef Dune::Fem::LevelGridPart< GridType > LevelGridPartType;
    LevelGridPartType levelGridPart( grid, level );
    typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, LevelGridPartType, 1 > LevelDiscreteFunctionSpaceType;
    LevelDiscreteFunctionSpaceType levelSpace( levelGridPart );
    typedef Dune::Fem::AdaptiveDiscreteFunction< LevelDiscreteFunctionSpaceType > LevelFunctionType;
    LevelFunctionType approximation( "approximation", levelSpace );
    interpolate( gridFunctionAdapter( Dune::Fem::ExactSolution<FunctionSpaceType>(), levelGridPart, levelSpace.order() + 2 ), approximation );

    typedef ToGridPart< LevelFunctionType, GridPartType > ToReference;
    Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
    Dune::Fem::H1Norm< GridPartType > h1norm( gridPart );
    error.push_back( { l2norm.distance( ToReference(approximation), reference ),
                       h1norm.distance( reference, ToReference(approximation) ) } );
  }
  for( int step = 1; step < error.size()-1; ++step )
  {
    double l2eoc = log( error[ step ][ 0 ] / error[ step -1 ][ 0 ] ) / log( 0.5 );
    double h1eoc = log( error[ step ][ 1 ] / error[ step -1 ][ 1 ] ) / log( 0.5 );

    // std::cout<< "L2 Eoc: " << l2eoc << " ( " << error[step-1][0] << " " << error[step][0] << " ) " << std::endl;
    // std::cout<< "H1 Eoc: " << h1eoc << " ( " << error[step-1][1] << " " << error[step][1] << " ) " << std::endl;

    if( l2eoc < space.order()+1 - 0.2 )
      DUNE_THROW(Dune::InvalidStateException,"EOC check for comparing with reference solution failed");
    if( h1eoc < space.order() - 0.2 )
      DUNE_THROW(Dune::InvalidStateException,"EOC check for comparing with reference solution failed");
  }
  if( error.back()[0] > 1e-12 )
    DUNE_THROW(Dune::InvalidStateException,"Final error not zero!" );
  if( error.back()[1] > 1e-12 )
    DUNE_THROW(Dune::InvalidStateException,"Final error not zero!");
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
