#ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_L2PROJECTION_HH
#define DUNE_FEM_OPERATOR_PROJECTION_LOCAL_L2PROJECTION_HH

#include <cstddef>

#include <type_traits>
#include <utility>

#include <dune/common/dynvector.hh>
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalL2Projection
    // -----------------

    /** \class LocalL2Projection
     *
     *  \brief please doc me
     *
     *  \tparam  BasisFunctionSet  basis function set
     */
    template< class BasisFunctionSet, class Implementation >
    class LocalL2Projection
    {
   public:
      /** \brief basis function set type */
      typedef BasisFunctionSet BasisFunctionSetType;

    protected:
#ifndef DOXYGEN

      LocalL2Projection () {}

#endif // #ifndef DOXYGEN

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      LocalL2Projection ( const LocalL2Projection & ) = default;

      /** \brief move constructor */
      LocalL2Projection ( LocalL2Projection && ) {}

      /** \brief assignment operator */
      LocalL2Projection &operator= ( const LocalL2Projection & ) = default;

      /** \brief move assignment operator */
      LocalL2Projection &operator= ( LocalL2Projection && ) {}

      /** \} */

    public:
      /** \name Public member methods
       *  \{
       */

      /** \brief return basis function set */
      BasisFunctionSet basisFunctionSet () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().basisFunctionSet() );
        return impl().basisFunctionSet();
      }

      /** \brief please doc me
       *
       *  \tparam  LocalFunction  local function type
       *  \tparam  LocalDofVector  local dof vector type
       *
       *  \param[in]  localFunction  local function
       *  \param[out]  localDofVector  dof vector
       */
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        apply( localFunction, localDofVector );
      }

      /** \brief please doc me
       *
       *  \tparam  LocalFunction  local function type
       *  \tparam  LocalDofVector  local dof vector type
       *
       *  \param[in]  localFunction  local function
       *  \param[out]  localDofVector  dof vector
       */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().apply( localFunction, localDofVector ) );
        impl().apply( localFunction, localDofVector );
      }

      /** \} */

    protected:
      const Implementation &impl () const
      {
        return static_cast< const Implementation & >( *this );
      }
    };



    // DefaultLocalL2Projection
    // ------------------------

    template< class LocalRieszProjection, class Quadrature >
    class DefaultLocalL2Projection
    : public LocalL2Projection< typename LocalRieszProjection::BasisFunctionSetType, DefaultLocalL2Projection< LocalRieszProjection, Quadrature > >
    {
      typedef DefaultLocalL2Projection< LocalRieszProjection, Quadrature > ThisType;
      typedef LocalL2Projection< typename LocalRieszProjection::BasisFunctionSetType, DefaultLocalL2Projection< LocalRieszProjection, Quadrature > > BaseType;

    public:
      /** \brief type of local Riesz type */
      typedef LocalRieszProjection LocalRieszProjectionType;
      /** \copydoc Dune::Fem::LocalL2Projection::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

    private:
      typedef typename BasisFunctionSetType::RangeType RangeType;
      typedef typename RangeType::value_type RangeFieldType;

    public:
      /** \name Construction
       *  \{
       */

      template< class... Args >
      explicit DefaultLocalL2Projection ( Args &&... args )
        : rieszProjection_( std::forward< Args >( args )... )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      DefaultLocalL2Projection ( const ThisType & ) = default;

#ifndef DOXYGEN

      DefaultLocalL2Projection ( ThisType &other )
        : DefaultLocalL2Projection( static_cast< const ThisType & >( other ) )
      {}

#endif // #ifndef DOXYGEN

      DefaultLocalL2Projection ( ThisType &&other )
        : rieszProjection_( std::move( other.rieszProjection_ ) )
      {}

      DefaultLocalL2Projection &operator= ( const ThisType & ) = default;

      DefaultLocalL2Projection &operator= ( ThisType &&other )
      {
        rieszProjection_ = std::move( other.rieszProjection_ );
        return *this;
      }

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** unbind */
      void unbind() {}

      /** \copydoc Dune::Fem::LocalL2Projection::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet () const
      {
        return rieszProjection_.basisFunctionSet();
      }

      /** \copydoc Dune::Fem::LocalL2Projection::apply */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        static_assert( std::is_same< RangeType, typename LocalFunction::RangeType >::value,
                       "RangeType and LocalFunction::RangeType have to be the same type" );

        const BasisFunctionSetType basisFunctionSet = this->basisFunctionSet();
        const auto geometry = basisFunctionSet.entity().geometry();

        f_.resize( basisFunctionSet.size() );
        f_ = static_cast< RangeFieldType >( 0 );

        Quadrature quadrature( geometry.type(), localFunction.order() + basisFunctionSet.order() );
        const std::size_t nop = quadrature.nop();
        for( std::size_t qp = 0; qp < nop; ++qp )
        {
          RangeType value;
          localFunction.evaluate( quadrature[ qp ], value );
          value *= quadrature.weight( qp )*geometry.integrationElement( quadrature.point( qp ) );
          basisFunctionSet.axpy( quadrature[ qp ], value, f_ );
        }

        rieszProjection_.apply( f_, dofs );
      }

    private:
      LocalRieszProjectionType rieszProjection_;
      mutable Dune::DynamicVector< RangeFieldType > f_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_L2PROJECTION_HH
