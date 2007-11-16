#ifndef DUNE_FEM_TEST_DGL2PROJECTION_HH
#define DUNE_FEM_TEST_DGL2PROJECTION_HH

#include <dune/fem/quadrature/cachequad.hh>

namespace Dune
{

  template< class DiscreteFunctionImp >
  class DGL2Projection
  {
  public:
    typedef DiscreteFunctionImp DiscreteFunctionType;

  private:
    typedef DGL2Projection< DiscreteFunctionType > ThisType;

  public:
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;
    
  public:
    template< class FunctionType >
    static void project( const FunctionType &function,
                         DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
      
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      typedef typename IteratorType :: Entity EntityType;

      typedef typename EntityType :: Geometry GeometryType;

      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
      
      const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();

      discreteFunction.clear();

      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;

        const GeometryType &geometry = entity.geometry();

        LocalFunctionType localFunction = discreteFunction.localFunction( entity );
        const unsigned int numDofs = localFunction.numDofs();
        const BaseFunctionSetType &baseFunctionSet = localFunction.baseFunctionSet();

        QuadratureType quadrature( entity, 2*dfSpace.order() );
        const unsigned int numQuadraturePoints = quadrature.nop();
        for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          RangeType phi;
          function.evaluate( geometry.global( quadrature.point( qp ) ), phi );

          const RangeFieldType weight = quadrature.weight( qp );

          for( unsigned int i = 0; i < numDofs; ++i )
          {
            RangeType psi;
            baseFunctionSet.evaluate( i, quadrature[ qp ], psi );
            localFunction[ i ] += weight * (phi * psi);
          }
        }
      }
    }
  };

}

#endif
