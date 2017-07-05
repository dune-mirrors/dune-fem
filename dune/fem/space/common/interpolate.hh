#ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
#define DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH

#include <vector>
#include <type_traits>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/space/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    // interpolate
    // -----------

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

    template< class Function, class DiscreteFunction, unsigned int partitions >
    static inline std::enable_if_t< !std::is_convertible< Function, HasLocalFunction >::value >
    interpolate ( const Function &u, DiscreteFunction &v, PartitionSet< partitions > ps )
    {
      typedef typename DiscreteFunction :: DiscreteFunctionSpaceType :: GridPartType  GridPartType;
      typedef GridFunctionAdapter< Function, GridPartType > GridFunctionType;

      GridFunctionType uGrid( "uGrid", u, v.space().gridPart() );

      interpolate( uGrid, v, ps );
    }

    template< class GridFunction, class DiscreteFunction, unsigned int partitions >
    static inline std::enable_if_t< std::is_convertible< GridFunction, HasLocalFunction >::value && Capabilities::hasInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType >::v >
    interpolate ( const GridFunction &u, DiscreteFunction &v, PartitionSet< partitions > ps )
    {
      ConstLocalFunction< GridFunction > uLocal( u );
      LocalContribution< DiscreteFunction, Assembly::Set > vLocal( v );

      // iterate over selected partition
      for( const auto entity : elements( v.gridPart(), ps ) )
      {
        // initialize u to entity
        uLocal.init( entity );

        // bind v to entity
        auto vGuard = bindGuard( vLocal, entity );

        // perform local interpolation
        v.space().interpolation( entity )( uLocal, vLocal );
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
