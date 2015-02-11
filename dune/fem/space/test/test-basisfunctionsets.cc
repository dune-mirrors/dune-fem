#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/fem/space/shapefunctionset/lagrange.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/simple.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>

#include "checkbasisfunctionset.hh"


// polynom approximation order of quadratures,
// at least polynom order of basis functions
const int polOrd = POLORDER;

//! the index set we are using
typedef Dune::GridSelector::GridType MyGridType;
typedef Dune::Fem::LeafGridPart< MyGridType > GridPartType;

static const int dimworld = MyGridType::dimensionworld;

// see dune/common/functionspace.hh
typedef Dune::Fem::FunctionSpace< MyGridType::ctype, double, dimworld, 1 > ScalarFunctionSpaceType;
typedef Dune::Fem::FunctionSpace< MyGridType::ctype, double, dimworld, dimworld > VectorialFunctionSpaceType;

static const int dimDomain = ScalarFunctionSpaceType::dimDomain;
static const int dimRange = 1;

void checkBasisFunctionSets ( const GridPartType &gridPart )
{
  typedef GridPartType::template Codim< 0 >::EntityType EntityType;

  // we only check on one element
  auto iterator = gridPart.template begin< 0 >();
  iterator++;
  const EntityType &entity = *iterator;

  // create quadrature
  const Dune::Fem::CachingQuadrature< GridPartType, 0 > quadrature( entity, polOrd );

  // needs a geometry type to construct
  typedef Dune::Fem::LagrangeShapeFunctionSet< ScalarFunctionSpaceType, polOrd > ScalarLagrangeShapeFunctionSetType;

  // needs an order to construct
  typedef Dune::Fem::LegendreShapeFunctionSet< ScalarFunctionSpaceType > ScalarLegendreShapeFunctionSetType;

  // needs a geometry type to construct
  typedef Dune::Fem::OrthonormalShapeFunctionSet< ScalarFunctionSpaceType, polOrd > ScalarOrthonormalShapeFunctionSetType;

  // type of error
  typedef Dune::FieldVector< double, 5 > ErrorType;

  // prepare shapefunctions
  ScalarLagrangeShapeFunctionSetType scalarLagrangeShapeFunctionSet( entity.type() );
  ScalarLegendreShapeFunctionSetType scalarLegendreShapeFunctionSet( polOrd );
  ScalarOrthonormalShapeFunctionSetType scalarOrthonormalShapeFunctionSet( entity.type() );

  ErrorType error;
  // default basis function set
  {
    Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLagrangeShapeFunctionSetType >
    basisSet1( entity, scalarLagrangeShapeFunctionSet );
    error = Dune::Fem::checkQuadratureConsistency( basisSet1, quadrature, true );
    if( error.two_norm() > 1e-8 )
      DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LagrangeShapeFunctionSet > test failed." );

    Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLegendreShapeFunctionSetType >
    basisSet2( entity, scalarLegendreShapeFunctionSet );
    error = Dune::Fem::checkQuadratureConsistency( basisSet2, quadrature, true );
    if( error.two_norm() > 1e-8 )
      DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LegendreShapeFunctionSet > test failed." );

    Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarOrthonormalShapeFunctionSetType >
    basisSet3( entity, scalarOrthonormalShapeFunctionSet );
    error = Dune::Fem::checkQuadratureConsistency( basisSet3, quadrature, false );
    if( error.two_norm() > 1e-8 )
      DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< OrthonormalShapeFunctionSet > test failed." );
  }

  // simple basis function set
  {}

  // vectorial basis function set
  {
    typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLegendreShapeFunctionSetType > ScalarBasisFunctionSetType;
    ScalarBasisFunctionSetType scalarBasisSet( entity, scalarLegendreShapeFunctionSet );

    Dune::Fem::VectorialBasisFunctionSet< ScalarBasisFunctionSetType,
                                          typename VectorialFunctionSpaceType::RangeType, Dune::Fem::VerticalDofAlignment >
    basisSet1( scalarBasisSet );
    error = Dune::Fem::checkQuadratureConsistency( basisSet1, quadrature, true );
    if( error.two_norm() > 1e-8 )
      DUNE_THROW( Dune::InvalidStateException, " VectorialBasisFunctionSet< ShapeFunctionSet, VerticalDofAlignment > test failed." );

    Dune::Fem::VectorialBasisFunctionSet< ScalarBasisFunctionSetType,
                                          typename VectorialFunctionSpaceType::RangeType, Dune::Fem::HorizontalDofAlignment >
    basisSet2( scalarBasisSet );
    error = Dune::Fem::checkQuadratureConsistency( basisSet2, quadrature, true );
    if( error.two_norm() > 1e-8 )
      DUNE_THROW( Dune::InvalidStateException, " VectorialBasisFunctionSet< ShapeFunctionSet, HorizontalDofAlignment > test failed." );
  }
}


//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  try
  {
    std::stringstream tmp;
    tmp << dimworld;
    std::string macroGridName( tmp.str());
    macroGridName += "dgrid.dgf";

    Dune::GridPtr< MyGridType > gridptr( macroGridName );
    MyGridType &grid = *gridptr;

    grid.globalRefine( 1 );

    GridPartType part( grid );

    if( part.template begin< 0 >() == part.template end< 0 >() )
      DUNE_THROW( Dune::InvalidStateException, "Macro grid has zero elements!" );

    checkBasisFunctionSets( part );

    return 0;
  }
  catch( const Dune::Exception &exception )
  {
    std::cerr << exception << std::endl;
    return 1;
  }
}

