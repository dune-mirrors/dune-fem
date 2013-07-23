// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef MASSOPERATOR_HH
#define MASSOPERATOR_HH

#include <config.h>

#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/misc/l2norm.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/file/datawriter.hh>

#if defined HAVE_PETSC

#include <dune/fem/operator/linear/petscoperator.hh>

class VariableFilenameParameter
: public Dune::Fem::DataOutputParameters 
{

public:
  VariableFilenameParameter ( const std::string &prefix ) 
  : prefix_( prefix )
  {}

  virtual std::string prefix () const
  {
    return prefix_;
  }
  
private:
  std::string prefix_;
};


template< class DiscreteFunction >
class MassOperator
: public Dune::Operator< typename DiscreteFunction::RangeFieldType,
                         typename DiscreteFunction::RangeFieldType,
                         DiscreteFunction, DiscreteFunction >
{
  typedef MassOperator< DiscreteFunction > ThisType;

public:
  typedef DiscreteFunction DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef Dune::Fem::DofManager< typename DiscreteFunctionSpaceType::GridType > DofManagerType;
  
  typedef Dune::Fem::CachingQuadrature< typename DiscreteFunctionSpaceType::GridPartType, 0 >
    QuadratureType;

  #if PETSCLINEAROPERATOR == 1
    typedef Dune::Fem::PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType >
      LinearOperatorType;
    #warning using Petsc for linear operator
  #else
    typedef Dune::Fem::SparseRowLinearOperator
      < DiscreteFunctionType, DiscreteFunctionType >
      LinearOperatorType;
    #warning using SparseRow for linear operator
  #endif


  explicit MassOperator ( const DiscreteFunctionSpaceType &dfSpace )
  : dfSpace_( dfSpace ),
    dofManager_( DofManagerType::instance( dfSpace.gridPart().grid() ) ),
    systemMatrix_( "mass operator", dfSpace, dfSpace ),
    sequence_( -1 )
  {}

  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
  {
    systemMatrix().apply( u, w );
  }

  template< class Function >
  void assembleRHS( const Function &u, DiscreteFunctionType &w ) const;

  const LinearOperatorType &systemMatrix () const;

private:
  const DiscreteFunctionSpaceType &dfSpace_;
  const DofManagerType &dofManager_;
  mutable LinearOperatorType systemMatrix_;
  mutable int sequence_;
};


template< class DiscreteFunction >
template< class Function >
void MassOperator< DiscreteFunction >
  ::assembleRHS( const Function &u, DiscreteFunctionType &w ) const
{
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType FieldType;
  typedef typename DiscreteFunctionSpaceType::RangeType      RangeType;

  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef ::Dune::Fem::L2Norm< typename DiscreteFunctionType::GridPartType > L2NormType;

  typedef Dune::tuple< DiscreteFunction* > IOTupleType;
  typedef Dune::Fem::DataOutput< Dune::GridSelector::GridType, IOTupleType > DataOutputType;

  w.clear();

  // run over entities
  const IteratorType end = dfSpace_.end();
  for( IteratorType it = dfSpace_.begin(); it != end; ++it )
  {

    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    LocalFunctionType localFunction = w.localFunction( entity );

    // run over quadrature points
    QuadratureType quadrature( entity, 2*dfSpace_.order()+1 );
    const unsigned int numQuadraturePoints = quadrature.nop();
    for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
    {
      // evaluate u
      const typename QuadratureType::CoordinateType &x = quadrature.point( qp );        

      RangeType uValue;
      u.evaluate( geometry.global( x ), uValue );

      // put all things together and don't forget quadrature weights
      const FieldType weight = quadrature.weight( qp )*geometry.integrationElement( x );

      // apply weight 
      uValue *= weight; 

      // add to local function 
      localFunction.axpy( quadrature[ qp ], uValue );
    }

  }

  w.communicate();

}


template< class DiscreteFunction >
const typename MassOperator< DiscreteFunction >::LinearOperatorType &
MassOperator< DiscreteFunction >::systemMatrix () const
{
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType FieldType;

  typedef typename LinearOperatorType::LocalMatrixType LocalMatrixType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  if( sequence_ != dofManager_.sequence() )
  {
    systemMatrix_.reserve(Dune::Fem::SimpleStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>());
    systemMatrix_.clear();

    std::vector< typename DiscreteFunctionSpaceType::RangeType > values;

    // run over entities
    const IteratorType end = dfSpace_.end();
    for( IteratorType it = dfSpace_.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      LocalMatrixType localMatrix = systemMatrix_.localMatrix( entity, entity );

      const BasisFunctionSetType &basis = localMatrix.domainBasisFunctionSet();
      const unsigned int numBasisFunctions = basis.size();
      values.resize( numBasisFunctions );

      // run over quadrature points
      QuadratureType quadrature( entity, 2*dfSpace_.order() );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        // evaluate base functions
        basis.evaluateAll( quadrature[ qp ], values );

        // get quadrature weight
        const typename QuadratureType::CoordinateType &x = quadrature.point( qp );
        const FieldType weight = quadrature.weight( qp )
                                    * geometry.integrationElement( x );

        // update system matrix
        for( unsigned int i = 0; i < numBasisFunctions; ++i )
        {
          // add diagonal entries
          localMatrix.add( i, i, weight * (values[ i ] * values[ i ]) );

          // add non-diagonal entries, use symmetric structure of the matrix
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
#if PETSCLINEAROPERATOR == 1
  systemMatrix_.communicate();
#endif

  return systemMatrix_;
}

#endif // #if defined HAVE_PETSC

#endif // #ifndef MASSOPERATOR_HH
