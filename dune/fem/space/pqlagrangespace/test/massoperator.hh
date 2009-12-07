#ifndef MASSOPERATOR_HH
#define MASSOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune
{

  template< class DiscreteFunction >
  class MassOperator
  : public Dune::Operator< typename DiscreteFunction::RangeFieldType,
                           typename DiscreteFunction::RangeFieldType,
                           DiscreteFunction, DiscreteFunction >
  {
    typedef MassOperator< DiscreteFunction > ThisType;

    struct MatrixTraits;

  public:
    typedef DiscreteFunction                                         DiscreteFunctionType;
    typedef typename DiscreteFunctionType
              :: DiscreteFunctionSpaceType                           DiscreteFunctionSpaceType;
    typedef DofManager< typename DiscreteFunctionSpaceType
                       :: GridType >                                 DofManagerType;
    
    typedef Dune :: CachingQuadrature
                      < typename DiscreteFunctionSpaceType
                          :: GridPartType, 0 >                       QuadratureType;

    typedef Dune :: SparseRowMatrixObject< DiscreteFunctionSpaceType,
                                           DiscreteFunctionSpaceType,
                                           MatrixTraits >            SystemMatrixType;

    explicit MassOperator ( const DiscreteFunctionSpaceType &dfSpace )
    : dfSpace_( dfSpace ),
      dofManager_( DofManagerType::instance( dfSpace.grid() ) ),
      systemMatrix_( dfSpace, dfSpace ),
      sequence_( -1 )
    {}

    virtual void
    operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
    {
      systemMatrix().apply( u, w );
    }

    template< class Function >
    void operator() ( const Function &u, DiscreteFunctionType &w ) const;

    const SystemMatrixType &systemMatrix () const;

  private:
    const DiscreteFunctionSpaceType &dfSpace_;
    const DofManagerType            &dofManager_;
    mutable SystemMatrixType         systemMatrix_;
    mutable int                      sequence_;
  };


  template< class DiscreteFunction >
  struct MassOperator< DiscreteFunction >::MatrixTraits
  {
    typedef DiscreteFunctionSpaceType                                RowSpaceType;
    typedef DiscreteFunctionSpaceType                                ColumnSpaceType;

    typedef LagrangeMatrixSetup< false >                             StencilType;

    typedef ParallelScalarProduct< DiscreteFunctionSpaceType >       ParallelScalarProductType;
  };


  template< class DiscreteFunction >
  template< class Function >
  void MassOperator< DiscreteFunction >
    ::operator() ( const Function &u, DiscreteFunctionType &w ) const
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType       IteratorType;
    typedef typename DiscreteFunctionSpaceType
              :: BaseFunctionSetType                                 BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType     FieldType;

    typedef typename DiscreteFunctionType :: LocalFunctionType       LocalFunctionType;
    typedef typename IteratorType :: Entity                          EntityType;
    typedef typename EntityType :: Geometry                          GeometryType;

    w.clear();
    std :: vector< typename DiscreteFunctionSpaceType :: RangeType > values;

    const IteratorType end = dfSpace_.end();
    for( IteratorType it = dfSpace_.begin(); it != end; ++it )
    {
      const EntityType &entity     = *it;
      const GeometryType &geometry = entity.geometry();

      LocalFunctionType localFunction = w.localFunction( entity );

      const BaseFunctionSetType &basis     = localFunction.baseFunctionSet();
      const unsigned int numBasisFunctions = basis.numBaseFunctions();
      values.resize( numBasisFunctions );

      QuadratureType quadrature( entity, 2*dfSpace_.order()+1 );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        const typename QuadratureType::CoordinateType &x = quadrature.point( qp );
        for( unsigned int i = 0; i < numBasisFunctions; ++i )
          basis.evaluate( i, quadrature[ qp ], values[ i ] );
        typename DiscreteFunctionSpaceType::RangeType uValue;
        u.evaluate( geometry.global( x ), uValue );

        const FieldType weight = quadrature.weight( qp )*geometry.integrationElement( x );
        for( unsigned int i = 0; i < numBasisFunctions; ++i )
          localFunction[ i ] += weight * (uValue*values[ i ]);
      }
    }
  }


  template< class DiscreteFunction >
  const typename MassOperator< DiscreteFunction >::SystemMatrixType &
  MassOperator< DiscreteFunction >::systemMatrix () const
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType       IteratorType;
    typedef typename DiscreteFunctionSpaceType
              :: BaseFunctionSetType                                 BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType     FieldType;

    typedef typename SystemMatrixType :: LocalMatrixType             LocalMatrixType;
    typedef typename IteratorType :: Entity                          EntityType;
    typedef typename EntityType :: Geometry                          GeometryType;

    if( sequence_ != dofManager_.sequence() )
    {
      systemMatrix_.reserve();
      systemMatrix_.clear();

      std::vector< typename DiscreteFunctionSpaceType::RangeType > values;

      const IteratorType end = dfSpace_.end();
      for( IteratorType it = dfSpace_.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        const GeometryType &geometry = entity.geometry();

        LocalMatrixType localMatrix = systemMatrix_.localMatrix( entity, entity );

        const BaseFunctionSetType &basis     = localMatrix.domainBaseFunctionSet();
        const unsigned int numBasisFunctions = basis.numBaseFunctions();
        values.resize( numBasisFunctions );

        QuadratureType quadrature( entity, 2*dfSpace_.order() );
        const unsigned int numQuadraturePoints = quadrature.nop();
        for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          const typename QuadratureType::CoordinateType &x = quadrature.point( qp );
          for( unsigned int i = 0; i < numBasisFunctions; ++i )
            basis.evaluate( i, quadrature[ qp ], values[ i ] );

          const FieldType weight = quadrature.weight( qp )*geometry.integrationElement( x );
          for( unsigned int i = 0; i < numBasisFunctions; ++i )
          {
            localMatrix.add( i, i, weight * (values[ i ] * values[ i ]) );
            for( unsigned int j = 0; j < i; ++j )
            {
              const FieldType value = weight * (values[ i ] * values[ j ]);
              localMatrix.add( i, j, value );
              localMatrix.add( j, i, value );
            }
          }
        }
      }
      sequence_ = dofManager_.sequence();
    }
    return systemMatrix_;
  }

}

#endif
