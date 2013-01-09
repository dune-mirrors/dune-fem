#ifndef DUNE_FEM_SPACE_INTERPOLATION_HH
#define DUNE_FEM_SPACE_INTERPOLATION_HH

namespace Dune
{

  namespace Fem
  {

    /** \class   Interpolation
     *  \ingroup DiscreteFunctionSpace
     *  \brief   native interpolation of a discrete function space
     *
     *  By definition of its degrees of freedom, each discrete function space
     *  has a native interpolation, which can be computed very quickly.
     *
     *  For example, the native interpolation of a Lagrange discrete function
     *  space is the evaluation in its Lagrange points.
     *  An orthonormal DG space would instead perform an \f$L^2\f$-Projection.
     *
     *  The actual implementation must locally be provided by the discrete
     *  function space through the method
     *  \code
template< class LocalFunction, class LocalDofVector >
void interpolate ( const LocalFunction &f, LocalDofVector &dofs ) const;
     *  \endcode
     *
     *  \tparam  DiscreteFunction  type of discrete function to interpolate to
     */
    template< class DiscreteFunction >
    struct Interpolation
    {
      typedef DiscreteFunction DiscreteFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      /** \brief interpolate a grid function
       *
       *  \param[in]   u  grid function to interpolate
       *  \param[out]  w  discrete function to represent the interpolation
       */
      template< class GridFunction >
      void operator() ( const GridFunction &u, DiscreteFunctionType &w )
      {
        typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
        typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

        typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

        const DiscreteFunctionSpaceType &space = w.space();

        const IteratorType end = space.end();
        for( IteratorType it = space.begin(); it != end; ++it )
        {
          const EntityType &entity = *it;
          LocalFunctionType wLocal = w.localFunction( entity );
          space.interpolate( u.localFunction( entity ), wLocal );
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_INTERPOLATION_HH
