#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_SINGLEVALUED_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_SINGLEVALUED_HH

#include <cassert>

#include <dune/fem/function/localfunction/average.hh>

#include "continuous.hh"

namespace Dune
{

  namespace Fem
  {

    // External forward declaration
    // ----------------------------

    template< class LocalFunction, class GridPart >
    class LocalAverage;



    // SingleValuedLocalFunction
    // -------------------------

    template< class Entity, class Range >
    class SingleValuedLocalFunction
    : public DefaultContinuousLocalFunction< Entity, Range, SingleValuedLocalFunction< Entity, Range > >
    {
      using ThisType = SingleValuedLocalFunction< Entity, Range >;
      using BaseType = DefaultContinuousLocalFunction< Entity, Range, SingleValuedLocalFunction< Entity, Range > >;

    public:
      /** \copydoc Dune::Fem::ContinuousLocalFunction::DomainType */
      using DomainType = typename BaseType::DomainType;
      /** \copydoc Dune::Fem::ContinuousLocalFunction::RangeType */
      using RangeType = typename BaseType::RangeType;
      /** \copydoc Dune::Fem::ContinuousLocalFunction::JacobianRangeType */
      using JacobianRangeType = typename BaseType::JacobianRangeType;
      /** \copydoc Dune::Fem::ContinuousLocalFunction::HessianRangeType */
      using HessianRangeType = typename BaseType::HessianRangeType;

      /** \copydoc Dune::Fem::ContinuousLocalFunction::EntityType */
      using EntityType = typename BaseType::EntityType;

      /** \name Construction
       *  \{
       */

      SingleValuedLocalFunction ()
        : entity_( nullptr )
      {}

      SingleValuedLocalFunction ( const EntityType &entity, const RangeType &value )
        : entity_( &entity ),
          value_( value )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      SingleValuedLocalFunction ( const ThisType & ) = default;

      /** \brief move constructor */
      SingleValuedLocalFunction ( ThisType && ) = default;

      /** \brief assignment operator */
      ThisType &operator= ( const ThisType & ) = default;

      /** \brief move assignment operator */
      ThisType &operator= ( ThisType && ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::ContinuousLocalFunction::order */
      static constexpr int order () noexcept { return 0; }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::evaluate */
      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        value = value_;
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::jacobian */
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 0 );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::hessian */
      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::entity */
      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      /** \} */

      /** \name Non-interface methods
       *  \{
       */

      explicit operator bool () const noexcept { return entity_; }

      /** \} */

    private:
      template< class, class > friend class LocalAverage;

      const EntityType *entity_;
      RangeType value_;
    };



    // Template specialization of LocalAverage< SingleValuedLocalFunction >
    // --------------------------------------------------------------------

    template< class GridPart, class Entity, class Range >
    class LocalAverage< SingleValuedLocalFunction< Entity, Range >, GridPart >
    {
      typedef SingleValuedLocalFunction< Entity, Range > LocalFunctionType;

    public:
      static void apply ( const LocalFunctionType &localFunction, typename LocalFunctionType::RangeType &average )
      {
        average = localFunction.value_;
      }

      void operator() ( const LocalFunctionType &localFunction, typename LocalFunctionType::RangeType &average ) const
      {
        apply( localFunction, average );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_SINGLEVALUED_HH
