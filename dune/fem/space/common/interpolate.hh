#ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
#define DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH

#include <dune/common/dynvector.hh>

#include <dune/grid/common/partitionset.hh>

namespace Dune
{

  namespace Fem
  {

    // interpolate
    // -----------

    template< class GridFunction, class DiscreteFunction, unsigned int partitions >
    static inline void interpolate ( const GridFunction &u, DiscreteFunction &v, PartitionSet< partitions > pset )
    {
      // obtain type of partition to iterator over
      const PartitionIteratorType pitype = derive_partition_iterator_type< partitions >::value;

      // obtain types
      typedef typename DiscreteFunction::GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename DiscreteFunction::GridPartType::template Codim< 0 >::template Partition< pitype >::IteratorType IteratorType;

      // reserve memory for local dof vector
      Dune::DynamicVector< typename DiscreteFunction::RangeFieldType > ldv;
      ldv.reserve( v.space().maxNumLocalDofs() );

      // iterate over selected partition of grid part
      const IteratorType end = v.gridPart().template end< 0, pitype >();
      for( IteratorType it = v.gridPart().template begin< 0, pitype >(); it != end; ++it )
      {
        // obtain entity
        const EntityType &entity = *it;

        // obtain local interpolation
        const auto interpolation = v.space().interpolation( entity );

        // resize local dof vector
        ldv.resize( v.space().basisFunctionSet( entity ).size() );

        // interpolate u locally
        interpolation( u.localFunction( entity ), ldv );

        // write local dofs into v
        v.setLocalDofs( entity, ldv );
      }
    }

    /**
     * \function interpolate
     * \ingroup  DiscreteFunctionSpace
     * \brief    perform native interpolation of a discrete function space
     *
     * By definition of its degrees of freedom, each discrete function space
     * has a native interpolation, which can be computed very quickly.
     *
     * For example, the native interpolation of a Lagrange discrete function
     * space is the evaluation in its Lagrange points.
     * An orthonormal DG space would instead perform an \f$L^2\f$-Projection.
     *
     * The actual implementation must locally be provided by the discrete
     * function space through the method
     * \code
     * template< class LocalFunction, class LocalDofVector >
     * void interpolate ( const LocalFunction &f, LocalDofVector &dofs ) const;
     * \endcode
     *
     * \param[in]   u  grid function to interpolate
     * \param[out]  v  discrete function to represent the interpolation
     */
    template< class GridFunction, class DiscreteFunction >
    static inline void interpolate ( const GridFunction &u, DiscreteFunction &v )
    {
      // just call interpolate for the all partition
      interpolate( u, v, Partitions::all );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
