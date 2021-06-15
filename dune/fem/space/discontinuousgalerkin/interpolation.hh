#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_INTERPOLATION_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_INTERPOLATION_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/agglomerationquadrature.hh>

#include <dune/fem/operator/projection/local/l2projection.hh>
#include <dune/fem/operator/projection/local/riesz.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class GridPart, class BasisFunctionSet,
              class Quadrature = CachingQuadrature< GridPart, BasisFunctionSet::EntityType::codimension > >
    class DiscontinuousGalerkinLocalL2Projection;

    template< class GridPart, class BasisFunctionSet,
              class Quadrature = CachingQuadrature< GridPart, BasisFunctionSet::EntityType::codimension >  >
    class LocalOrthonormalL2Projection;


    // LocalOrthonormalL2Projection
    // ----------------------------

    /** \brief specialization of local L2 projection for orthonormal DG spaces */
    template< class GridPart, class BasisFunctionSet, class Quadrature >
    class LocalOrthonormalL2Projection
    : public LocalL2Projection< BasisFunctionSet, LocalOrthonormalL2Projection< GridPart, BasisFunctionSet, Quadrature > >
    {
      typedef LocalOrthonormalL2Projection< GridPart, BasisFunctionSet, Quadrature > ThisType;
      typedef LocalL2Projection< BasisFunctionSet, LocalOrthonormalL2Projection< GridPart, BasisFunctionSet, Quadrature > > BaseType;

    public:
      /** \copydoc Dune::Fem::LocalL2Projection::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType      BasisFunctionSetType;
      /** \copydoc Dune::Fem::LocalL2Projection::EntityType */
      typedef typename BasisFunctionSetType :: EntityType  EntityType;

    private:
      typedef GridPart GridPartType;
      typedef typename GridPartType::GridType GridType;
      typedef typename BasisFunctionSetType :: RangeType RangeType;

    public:
      /** \name Construction
       *  \{
       */

      explicit LocalOrthonormalL2Projection ( const BasisFunctionSetType &basisFunctionSet )
        : basisFunctionSet_( basisFunctionSet )
      {}

      explicit LocalOrthonormalL2Projection ( BasisFunctionSetType &&basisFunctionSet )
        : basisFunctionSet_( std::forward< BasisFunctionSetType >( basisFunctionSet ) )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      LocalOrthonormalL2Projection ( const ThisType & ) = default;

      LocalOrthonormalL2Projection ( ThisType &&other ) = default;

      ThisType &operator= ( const ThisType & ) = default;

      ThisType &operator= ( ThisType &&other ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::LocalL2Projection::basisFunctionSet */
      const BasisFunctionSet& basisFunctionSet () const
      {
        return basisFunctionSet_;
      }

      /** \copydoc Dune::Fem::LocalL2Projection::apply */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        // get entity and geometry
        const EntityType &entity = localFunction.entity();

        if( entity.type().isNone() )
        {
          typedef ElementQuadrature< GridPartType, EntityType::codimension > ElementQuadratureType;
          ElementQuadratureType quadrature( entity, localFunction.order() + basisFunctionSet().order() );
          computeL2Projection( entity, quadrature, localFunction, localDofVector );
        }
        else
        {
          // create quadrature with appropriate order
          Quadrature quadrature( entity, localFunction.order() + basisFunctionSet().order() );
          computeL2Projection( entity, quadrature, localFunction, localDofVector );
        }
      }

      /** \} */

    protected:
      template <class QuadImpl, class LocalFunction, class LocalDofVector >
      void computeL2Projection( const EntityType& entity,
                                const QuadImpl& quadrature,
                                const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        // set all dofs to zero
        localDofVector.clear();

        const int nop = quadrature.nop();
        // adjust size of values
        values_.resize( nop );

        // evaluate local function for all quadrature points
        localFunction.evaluateQuadrature( quadrature, values_ );

        // apply weight only (for orthonormal basis set integration element and
        // mass matrix can be ignored even if geometry is non-affine)
        for(auto qp : quadrature )
          values_[ qp.index() ] *= qp.weight();

        // add values to local dof vector
        basisFunctionSet().axpy( quadrature, values_, localDofVector );
      }

      BasisFunctionSetType basisFunctionSet_;
      mutable std::vector< RangeType > values_;
    };


    // DiscontinuousGalerkinLocalL2Projection
    // --------------------------------------

    template< class GridPart, class BasisFunctionSet, class Quadrature >
    class DiscontinuousGalerkinLocalL2Projection
    : public LocalL2Projection< BasisFunctionSet, DiscontinuousGalerkinLocalL2Projection< GridPart, BasisFunctionSet, Quadrature > >
    {
      typedef DiscontinuousGalerkinLocalL2Projection< GridPart, BasisFunctionSet, Quadrature > ThisType;
      typedef LocalL2Projection< BasisFunctionSet, DiscontinuousGalerkinLocalL2Projection< GridPart, BasisFunctionSet, Quadrature > > BaseType;

    public:
      /** \copydoc Dune::Fem::LocalL2Projection::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

    private:
      typedef GridPart GridPartType;
      typedef typename GridPartType::GridType GridType;
      static const bool cartesian = Dune::Capabilities::isCartesian< GridType >::v;

      typedef typename std::conditional< cartesian,
          OrthonormalLocalRieszProjection< BasisFunctionSetType >,
          DenseLocalRieszProjection< BasisFunctionSetType, Quadrature >
        >::type LocalRieszProjectionType;

      typedef DefaultLocalL2Projection< LocalRieszProjectionType, Quadrature > Implementation;

    public:
      /** \name Construction
       *  \{
       */

      explicit DiscontinuousGalerkinLocalL2Projection ( const BasisFunctionSetType &basisFunctionSet )
        : impl_( LocalRieszProjectionType( basisFunctionSet ) )
      {}

      explicit DiscontinuousGalerkinLocalL2Projection ( BasisFunctionSetType &&basisFunctionSet )
        : impl_( LocalRieszProjectionType( std::forward< BasisFunctionSetType >( basisFunctionSet ) ) )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      DiscontinuousGalerkinLocalL2Projection ( const ThisType & ) = default;

      DiscontinuousGalerkinLocalL2Projection ( ThisType &&other )
        : impl_( std::move( other.impl_ ) )
      {}

      ThisType &operator= ( const ThisType & ) = default;

      ThisType &operator= ( ThisType &&other )
      {
        impl_ = std::move( other.impl_ );
        return *this;
      }

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::LocalL2Projection::basisFunctionSet */
      BasisFunctionSet basisFunctionSet () const
      {
        return impl_.basisFunctionSet();
      }

      /** \copydoc Dune::Fem::LocalL2Projection::apply */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        impl_( localFunction, localDofVector );
      }

      void unbind() {}

      /** \} */

    private:
      Implementation impl_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_INTERPOLATION_HH
