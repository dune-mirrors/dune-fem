#ifndef DUNE_FEM_SPACE_COMMON_LOCALINTERPOLATION_HH
#define DUNE_FEM_SPACE_COMMON_LOCALINTERPOLATION_HH

#include <optional>

#include <dune/fem/common/bindguard.hh>

namespace Dune
{

  namespace Fem
  {

    // Forward Declarations
    // --------------------

    template< class DiscreteFunctionSpace >
    class LocalInterpolation;

    template < class DiscreteFunctionSpace >
    class LocalInterpolation
    {
      // do not copy this class
      LocalInterpolation( const LocalInterpolation& ) = delete;
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    protected:
      typedef typename DiscreteFunctionSpaceType :: InterpolationType InterpolationType;

    public:
      typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;

      LocalInterpolation( const DiscreteFunctionSpaceType& space )
        : interpolation_( space.interpolation() )
      {}

      /** \brief initialize the local interpolation for an entity
       *
          Binds the local interpolation to an basisFunctionSet and entity.

          \param[in] entity to bind the local interpolation to
       */
      void bind ( const EntityType& entity )
      {
        interpolation_.bind( entity );
      }

      /** \brief clears the local interpolation by removing the basisFunctionSet
       *
       */
      void unbind ()
      {
        interpolation_.unbind();
      }

      /** \brief computes interpolation of locaFunction on entity and stores result in dofs

          \param[in]  localFunction object to be interpolated
          \param[out] dofs  vector to store result
       */
      template< class LocalFunction, class LocalDofVector >
      void operator () ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        interpolation_( localFunction, dofs );
      }

    protected:
      InterpolationType interpolation_;
    };


    template< class DiscreteFunctionSpace >
    class LocalInterpolationWrapper
    {
      typedef typename DiscreteFunctionSpace :: InterpolationImplType
        InterpolationImplType;

    public:
      typedef typename DiscreteFunctionSpace::EntityType  EntityType;

      LocalInterpolationWrapper( const DiscreteFunctionSpace& space )
        : space_( space )
        , interpolation_()
      {}

      void bind( const EntityType& entity )
      {
        interpolation_.emplace( space_.localInterpolation( entity ) );
      }

      void unbind()
      {
        interpolation().unbind();
      }

      /** \brief TODO, documentation */
      template< class LocalFunction, class LocalDofVector >
      void operator () ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        interpolation()( localFunction, dofs );
      }

    protected:
      const InterpolationImplType& interpolation() const
      {
        assert( interpolation_.has_value() );
        return *interpolation_;
      }
      InterpolationImplType& interpolation()
      {
        assert( interpolation_.has_value() );
        return *interpolation_;
      }

      const DiscreteFunctionSpace& space_;
      std::optional< InterpolationImplType > interpolation_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_LOCALINTERPOLATION_HH
