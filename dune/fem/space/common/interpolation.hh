#ifndef DUNE_FEM_SPACE_INTERPOLATION_HH
#define DUNE_FEM_SPACE_INTERPOLATION_HH

#include <dune/fem/misc/iteratorprovider.hh>

namespace Dune
{

  namespace Fem
  {

    // Interpolation
    // -------------

    /** \class   Interpolation
        \ingroup DiscreteFunctionSpace
        \brief   native interpolation of a discrete function space
      
        By definition of its degrees of freedom, each discrete function space
        has a native interpolation, which can be computed very quickly.
      
        For example, the native interpolation of a Lagrange discrete function
        space is the evaluation in its Lagrange points.
        An orthonormal DG space would instead perform an \f$L^2\f$-Projection.
      
        The actual implementation must locally be provided by the discrete
        function space through the method
        \code
template< class LocalFunction, class LocalDofVector >
void interpolate ( const LocalFunction &f, LocalDofVector &dofs ) const;
        \endcode
      
        \tparam  DiscreteFunction  type of discrete function to interpolate to
     */
    template< class DiscreteFunction,
              class IteratorProvider = Fem::IteratorProvider< typename DiscreteFunction::DiscreteFunctionSpaceType > >
    struct Interpolation
    {
      typedef DiscreteFunction DiscreteFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      /** \brief interpolate a grid function
        
          \param[in]   u  grid function to interpolate
          \param[out]  v  discrete function to represent the interpolation
       */
      template< class GridFunction >
      void operator() ( const GridFunction &u, DiscreteFunctionType &v )
      {
        apply( u, v );
      }

      /** \brief interpolate a grid function
        
          \param[in]   u  grid function to interpolate
          \param[out]  v  discrete function to represent the interpolation
       */
      template< class GridFunction >
      static void apply ( const GridFunction &u, DiscreteFunctionType &v )
      {
        typedef typename IteratorProvider::IteratorType IteratorType;
        typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
        typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

        const DiscreteFunctionSpaceType &space = v.space();
        IteratorProvider iteratorProvider( space );

        const IteratorType end = iteratorProvider.end();
        for( IteratorType it = iteratorProvider.begin(); it != end; ++it )
        {
          const EntityType &entity = *it;
          LocalFunctionType vLocal = v.localFunction( entity );
          space.interpolate( u.localFunction( entity ), vLocal );
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_INTERPOLATION_HH
