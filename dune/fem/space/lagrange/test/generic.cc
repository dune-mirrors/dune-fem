#include <config.h>

#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>

#include "../shapefunctionset.hh"
#include "../lagrangepoints.hh"

using namespace Dune;
using namespace Fem;

#define VERBOSITY_LEVEL 0

#ifndef DIMENSION
#define DIMENSION 2
#endif

#ifndef POLORDER
#define POLORDER 2
#endif

#ifndef TOPOLOGYTYPE
#define TOPOLOGYTYPE Dune::GeometeyTypes::simplex
#endif

typedef FunctionSpace< double, double, DIMENSION, 1 > FunctionSpaceType;

typedef FunctionSpaceType :: DomainType DomainType;
typedef FunctionSpaceType :: RangeType RangeType;

typedef Fem::LagrangeShapeFunctionFactory< FunctionSpaceType, POLORDER > ShapeFunctionFactoryType;
typedef ShapeFunctionFactoryType::ShapeFunctionType ShapeFunctionType;

typedef Fem::LagrangePoint< TOPOLOGYTYPE(DIMENSION).id(), DIMENSION, POLORDER > LagrangePointType;
typedef Fem::LagrangePointListImplementation< double, TOPOLOGYTYPE(DIMENSION).id(), DIMENSION, POLORDER >
  LagrangePointListType;


int main( int argc, char **argv )
{
  GeometryType geometryType = TOPOLOGYTYPE(DIMENSION);
  ShapeFunctionFactoryType shapeFunctionFactory ( geometryType );

  const unsigned int numShapeFunctions = shapeFunctionFactory.numShapeFunctions();

  std::cout << "Number of shape functions: " << numShapeFunctions;
  std::cout << std::endl << std::endl;

  Hybrid::forEach( std::make_index_sequence< DIMENSION+1 >{},
    [ & ]( auto i ){ std::cout << "MaxDofs< " << i << " >: " << LagrangePointType::template Codim< i >::maxDofs() << std::endl; } );
  std::cout << std::endl;

  unsigned int errors = 0;
  unsigned int indexErrors = 0;
  unsigned int pointSetErrors = 0;
  LagrangePointListType pointSet( 0 );
  for( unsigned int i = 0; i < numShapeFunctions; ++i )
  {
    LagrangePointType point( i );

    DomainType x;
    point.local( x );
    x -= pointSet.point( i );
    if( x.two_norm() >= 1e-14 )
      ++pointSetErrors;

    unsigned int codim, setCodim;
    unsigned int subEntity, setSubEntity;
    unsigned int dofNumber, setDofNumber;
    point.dofSubEntity( codim, subEntity, dofNumber );
    pointSet.dofSubEntity( i, setCodim, setSubEntity, setDofNumber );
    if( (setCodim != codim) || (setSubEntity != subEntity)
        || (setDofNumber != dofNumber) )
      ++pointSetErrors;

    unsigned int numDofs = point.numDofs( codim, subEntity );
    if( numDofs != pointSet.numDofs( codim, subEntity ) )
      ++pointSetErrors;

    unsigned int entityDofNumber
      = point.entityDofNumber( codim, subEntity, dofNumber );
    if( entityDofNumber != pointSet.entityDofNumber( codim, subEntity, dofNumber ) )
      ++pointSetErrors;

    if( dofNumber >= numDofs )
      ++errors;
    if( entityDofNumber != i )
      ++indexErrors;

    std :: cout << i << ": x = " << pointSet.point( i );
    std :: cout << ", height = " << point.height();
    std :: cout << ", codim = " << codim << ", subEntity = " << subEntity;
    std :: cout << ", dofNumber = " << dofNumber << " ( " << numDofs << " )";
    std :: cout << ", entityDofNumber = " << entityDofNumber;
    std :: cout << std :: endl;
  }
  std :: cout << std :: endl;
  std :: cout << "DoFs with dofNumber >= numDofs: " << errors << std :: endl;
  std :: cout << "DoFs entityDofNumber( dofSubEntity ) != i: " << indexErrors
              << std :: endl;
  std :: cout << "Errors in point set: " << pointSetErrors << std :: endl;
  std :: cout << std :: endl;

  errors = 0;
  for( unsigned int j = 0; j < numShapeFunctions; ++j )
  {
    ShapeFunctionType *shapeFunction = shapeFunctionFactory.createShapeFunction( j );

    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << "Shape function: " << j << std :: endl;
    #endif
    for( unsigned int i = 0; i < numShapeFunctions; ++i )
    {
      LagrangePointType point( i );

      DomainType x;
      point.local( x );

      RangeType phi;
      shapeFunction->evaluate( x, phi );

      double expected = ((i == j) ? 1.0 : 0.0);
      if( fabs( phi[ 0 ] - expected ) > 1e-8 )
        ++errors;

      #if (VERBOSITY_LEVEL >= 1)
        std :: cout << i << ": x = " << x;
        std :: cout << ", f( x ) = " << phi[ 0 ];
        std :: cout << std :: endl;
      #endif
    }

    delete shapeFunction;
    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << std :: endl;
    #endif
  }
  std :: cout << "Shape function evaluation errors: " << errors << std :: endl;

  errors = 0;
  for( unsigned int i = 0; i < numShapeFunctions; ++i )
  {
    LagrangePointType point( i );

    DomainType x;
    point.local( x );

    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << "Lagrange point " << i << ": " << x << std :: endl;
    #endif

    typedef ShapeFunctionType::JacobianRangeType JacobianRangeType;
    JacobianRangeType jacobian( 0 );
    for( unsigned int j = 0; j < numShapeFunctions; ++j )
    {
      ShapeFunctionType *shapeFunction = shapeFunctionFactory.createShapeFunction( j );

      JacobianRangeType tmp;
      shapeFunction->jacobian( x, tmp );
      #if (VERBOSITY_LEVEL >= 1)
        std :: cout << "ShapeFunction " << j << ", jacobian: " << tmp << std :: endl;
      #endif

      jacobian += tmp;

      delete shapeFunction;
    }

    if( fabs( jacobian[ 0 ].two_norm() ) > 1e-8 )
      ++errors;

    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << std :: endl;
    #endif
  }
  std :: cout << "Shape function derivtive summation errors: " << errors << std :: endl;

  return 0;
}
