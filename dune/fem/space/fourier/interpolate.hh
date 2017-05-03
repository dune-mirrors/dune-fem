#ifndef DUNE_FEM_SPACE_FOURIER_INTERPOLATE_HH
#define DUNE_FEM_SPACE_FOURIER_INTERPOLATE_HH

#include <limits>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/fourier/declaration.hh>
#include <dune/fem/space/fourier/functionset.hh>

namespace Dune
{

  namespace Fem
  {

    // interpolate
    // -----------

    template< class GridFunction, class DiscreteFunction, unsigned int partitions >
    static inline std::enable_if_t< std::is_convertible< GridFunction, HasLocalFunction >::value && IsFourierDiscreteFunctionSpace< typename DiscreteFunction::DiscreteFunctionSpaceType >::value >
    interpolate ( const GridFunction &u, DiscreteFunction &v, PartitionSet< partitions > ps )
    {
      typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

      const int dimDomain = DiscreteFunctionSpaceType::dimDomain;
      const int dimRange = DiscreteFunctionSpaceType::dimRange;

      const int size = FourierFunctionSetSize< dimDomain, DiscreteFunctionSpaceType::polynomialOrder-1 >::v;

      // mass matrix and right hand side for L2 projection
      FieldMatrix< RangeFieldType, size, size > mat( 0 );
      FieldMatrix< RangeFieldType, dimRange, size > rhs( 0 );

      const auto &functionSet = v.space().functionSet();
      if( functionSet.order() == std::numeric_limits< int >::max() )
        DUNE_THROW( InvalidStateException, "Cannot interpolate to Fourier space if quadrature order is not specified" );
      typename GridFunction::LocalFunctionType uLocal( u );

      // Note: The partition is ignored, here. As all degrees of freedom are
      //       global, we need to fill in all DoFs anyway.
      //       As we need to compute the global mass matrix, we only add our
      //       interior contibutions, though.
      for( const auto entity : elements( v.gridPart(), Partitions::interior ) )
      {
        const auto geometry = entity.geometry();

        uLocal.init( entity );

        for( const auto qp : Dune::Fem::CachingQuadrature< GridPartType, 0 >( entity, 2*functionSet.order() ) )
        {
          // obtain global position and weight
          const auto x = geometry.global( qp.position() );
          const RangeFieldType weight = geometry.integrationElement( x ) * qp.weight();

          // evaluate scalar function set
          FieldVector< RangeFieldType, size > values;
          functionSet.evaluateEach( x, [ &values ] ( std::size_t i, const RangeFieldType &v ) { values[ i ] = v; } );

          // update mass matrix
          for( int i = 0; i < size; ++i )
            mat[ i ].axpy( values[ i ]*weight, values );

          // evaluate u
          RangeType uValue;
          uLocal.evaluate( qp, uValue );

          // update right hand side
          for( int i = 0; i < dimRange; ++i )
            rhs[ i ].axpy( uValue[ i ]*weight, values );
        }
      }

      // globally sum up mat and rhs
      std::array< RangeFieldType, (size+dimRange)*size > buffer;
      for( int i = 0; i < size; ++i )
        std::copy_n( mat[ i ].begin(), size, buffer.begin() + i*size );
      for( int i = 0; i < dimRange; ++i )
        std::copy_n( rhs[ i ].begin(), size, buffer.begin() + (i+size)*size );
      v.gridPart().comm().sum( buffer.data(), buffer.size() );
      for( int i = 0; i < size; ++i )
        std::copy_n( buffer.begin() + i*size, size, mat[ i ].begin() );
      for( int i = 0; i < dimRange; ++i )
        std::copy_n( buffer.begin() + (i+size)*size, size, rhs[ i ].begin() );

      // solve mat * dofs^T = rhs^T
      mat.invert();
      FieldMatrix< RangeFieldType, dimRange, size > dofs( 0 );
      for( int i = 0; i < dimRange; ++i )
        mat.umv( rhs[ i ], dofs[ i ] );

      // copy dofs to dof vector of v
      for( int i = 0; i < size; ++i )
        for( int j = 0; j < dimRange; ++j )
          v.dofVector()[ 0 ][ i*dimRange + j ] = dofs[ j ][ i ];
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_INTERPOLATE_HH
