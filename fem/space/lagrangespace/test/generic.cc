#include <config.h>

#include "../basefunctions.hh"
#include "../lagrangepoints.hh"

using namespace Dune;

#define VERBOSITY_LEVEL 0

#ifndef DIMENSION
#define DIMENSION 2
#endif

#ifndef POLORDER
#define POLORDER 2
#endif

#ifndef GEOMETRYTYPE
#define GEOMETRYTYPE simplex
#endif

typedef FunctionSpace< double, double, DIMENSION, 1 > FunctionSpaceType;

typedef FunctionSpaceType :: DomainType DomainType;
typedef FunctionSpaceType :: RangeType RangeType;

typedef LagrangeBaseFunctionFactory< FunctionSpaceType, DIMENSION, POLORDER >
  BaseFunctionFactoryType;
typedef BaseFunctionInterface< FunctionSpaceType > BaseFunctionType;

typedef LagrangePoint< GeometryType :: GEOMETRYTYPE, DIMENSION, POLORDER >
  LagrangePointType;
typedef LagrangePointListImplementation
  < double, GeometryType :: GEOMETRYTYPE, DIMENSION, POLORDER >
  LagrangePointListType;



template< unsigned int codim >
class printMaxDofs
{
public:
  static void print()
  {
    const unsigned int maxCodimDofs
      = LagrangePointType :: template Codim< codim > :: maxDofs();
  
    std :: cout << "MaxDofs< " << codim << " >: " << maxCodimDofs << std :: endl;
    printMaxDofs< codim-1 > :: print();
  }
};



template<>
class printMaxDofs< 0 >
{
private:
  enum { codim = 0 };

public:
  static void print()
  {
    const unsigned int maxCodimDofs
      = LagrangePointType :: Codim< codim > :: maxDofs();
  
    std :: cout << "MaxDofs< " << codim << " >: " << maxCodimDofs << std :: endl;
  }
};



int main( int argc, char **argv )
{
  GeometryType geometryType( GeometryType :: GEOMETRYTYPE, DIMENSION );
  BaseFunctionFactoryType baseFunctionFactory ( geometryType );

  const unsigned int numBaseFunctions = baseFunctionFactory.numBaseFunctions();

  std :: cout << "Number of base functions: " << numBaseFunctions;
  std :: cout << std :: endl << std :: endl;

  printMaxDofs< DIMENSION > :: print();
  std :: cout << std :: endl;

  unsigned int errors = 0;
  unsigned int indexErrors = 0;
  unsigned int pointSetErrors = 0;
  LagrangePointListType pointSet( 0 );
  for( unsigned int i = 0; i < numBaseFunctions; ++i )
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
  for( unsigned int j = 0; j < numBaseFunctions; ++j )
  {
    BaseFunctionType &baseFunction = *(baseFunctionFactory.baseFunction( j ));

    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << "Base function: " << j << std :: endl;
    #endif
    for( unsigned int i = 0; i < numBaseFunctions; ++i )
    {
      LagrangePointType point( i );

      DomainType x;
      point.local( x );

      RangeType phi;
      FieldVector< deriType, 0 > derivative;
      baseFunction.evaluate( derivative, x, phi );
      
      double expected = ((i == j) ? 1.0 : 0.0);
      if( fabs( phi[ 0 ] - expected ) > 1e-8 )
        ++errors;

      #if (VERBOSITY_LEVEL >= 1)
        std :: cout << i << ": x = " << x;
        std :: cout << ", f( x ) = " << phi[ 0 ];
        std :: cout << std :: endl;
      #endif
    }

    delete &baseFunction;
    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << std :: endl;
    #endif
  }
  std :: cout << "Base function evaluation errors: " << errors << std :: endl;
  
  errors = 0;
  for( unsigned int i = 0; i < numBaseFunctions; ++i ) {
    LagrangePointType point( i );

    DomainType x;
    point.local( x );

    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << "Lagrange point " << i << ": " << x << std :: endl;
    #endif

    for( deriType k = 0; k < DIMENSION; ++k ) {
      FieldVector< deriType, 1 > derivative( k );

      RangeType sum( 0 );
      for( unsigned int j = 0; j < numBaseFunctions; ++j )
      {
        BaseFunctionType &baseFunction = *(baseFunctionFactory.baseFunction( j ));

        RangeType phi;
        baseFunction.evaluate( derivative, x, phi );
        #if (VERBOSITY_LEVEL >= 1)
          std :: cout << "BaseFunction " << j << ", derivative " << k << ": " << phi[ 0 ] << std :: endl;
        #endif
        sum += phi;

        delete &baseFunction;
      }

      #if (VERBOSITY_LEVEL >= 1)
        std :: cout << "Derivative " << k << " sums up to " << sum[ 0 ] << "." << std :: endl;
      #endif
      if( fabs( sum[ 0 ] ) > 1e-8 )
        ++errors;
    }
    #if (VERBOSITY_LEVEL >= 1)
      std :: cout << std :: endl;
    #endif
  }
  std :: cout << "Base function derivtive summation errors: " << errors << std :: endl;
  
  return 0;
}
