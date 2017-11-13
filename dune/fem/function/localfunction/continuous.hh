#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONTINUOUS_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_CONTINUOUS_HH

#include <cassert>
#include <cstddef>

#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // ContinuousLocalFunction
    // -----------------------

    template< class FunctionSpace, class Entity, class Implementation >
    class ContinuousLocalFunction
    {
    public:
      /** \brief function space type */
      using FunctionSpaceType = FunctionSpace;

      /** \brief domain field type */
      using DomainFieldType = typename FunctionSpaceType::DomainFieldType;
      /** \brief range field type */
      using RangeFieldType = typename FunctionSpaceType::RangeFieldType;

      /** \brief domain dimension */
      static const int dimDomain = FunctionSpaceType::dimDomain;
      /** \brief range dimension */
      static const int dimRange = FunctionSpaceType::dimRange;

      /** \brief domain type */
      using DomainType = typename FunctionSpaceType::DomainType;
      /** \brief range type */
      using RangeType = typename FunctionSpaceType::RangeType;
      /** \brief jacobian range type */
      using JacobianRangeType = typename FunctionSpaceType::JacobianRangeType;
      /** \brief hessian range type */
      using HessianRangeType = typename FunctionSpaceType::HessianRangeType;

      /** \brief entity type */
      using EntityType = Entity;
      /** \brief local coordinate type */
      using LocalCoordinateType = typename EntityType::Geometry::LocalCoordinate;

    protected:
      ContinuousLocalFunction () = default;

    public:
      /** \name Query methods
       *  \{
       */

      /** \brief obtain the order of this local function */
      int order () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().order() );
        return impl().order();
      }

      /** \} */

      /** \name Evaluation
       *  \{
       */

      /** \brief evaluate local function */
      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().evaluate( x, value ) );
        impl().evaluate( x, value );
      }

      /** \brief evaluate jacobian of local function */
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().jacobian( x, jacobian ) );
        impl().jacobian( x, jacobian );
      }

      /** \brief evaluate hessian of local function */
      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().hessian( x, hessian ) );
        impl().hessian( x, hessian );
      }

      /** \brief evaluate local function in given quadrature points */
      template< class Quadrature, class RangeArray >
      void evaluateQuadrature ( const Quadrature &quadrature, RangeArray &values ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().evaluateQuadrature ( quadrature, values ) );
        impl().evaluateQuadrature ( quadrature, values );
      }

      /** \} */

      /** \name Entity access
       *  \{
       */

      /** \brief obtain the entity, this local function lives on */
      const EntityType &entity () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().entity() );
        return impl().entity();
      }

      /** \} */

    protected:
      const Implementation &impl () const
      {
        return static_cast< const Implementation & >( *this );
      }
    };



    // DefaultContinuousLocalFunction
    // ------------------------------

    template< class Entity, class Range, class Implementation >
    class DefaultContinuousLocalFunction
    : public ContinuousLocalFunction< Dune::Fem::FunctionSpace< typename Entity::Geometry::ctype, typename Range::value_type, Entity::Geometry::coorddimension, Range::dimension >, Entity, Implementation >
    {
      using BaseType = ContinuousLocalFunction< Dune::Fem::FunctionSpace< typename Entity::Geometry::ctype, typename Range::value_type, Entity::Geometry::coorddimension, Range::dimension >, Entity, Implementation >;

    protected:
      using BaseType::impl;

    public:
      /** \brief range type */
      using RangeType = typename BaseType::RangeType;
      /** \brief jacobian range type */
      using JacobianRangeType = typename BaseType::JacobianRangeType;

      /** \name Evaluation
       *  \{
       */

      /** \copydoc Dune::Fem::ContinuousLocalFunction::evaluateQuadrature */
      template< class Quadrature, class RangeArray >
      void evaluateQuadrature ( const Quadrature &quadrature, RangeArray &values ) const
      {
        assert( quadrature.nop() > 0 );
        evaluateQuadrature( quadrature, values, &values[ 0 ] );
      }

      /** \} */

    protected:
      template< class Quadrature, class RangeArray >
      void evaluateQuadrature ( const Quadrature &quadrature, RangeArray &values, RangeType * ) const
      {
        const std::size_t nop = quadrature.nop();
        for( std::size_t qp = 0; qp < nop; ++qp )
          this->evaluate( quadrature[ qp ], values[ qp ] );
      }

      template< class Quadrature, class JacobianRangeArray >
      void evaluateQuadrature ( const Quadrature &quadrature, JacobianRangeArray &jacobians, JacobianRangeType * ) const
      {
        const std::size_t nop = quadrature.nop();
        for( std::size_t qp = 0; qp < nop; ++qp )
          this->jacobian( quadrature[ qp ], jacobians[ qp ] );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONTINUOUS_HH
