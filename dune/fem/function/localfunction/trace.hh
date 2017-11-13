#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_TRACE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_TRACE_HH

#include <cstddef>

#include <functional>

#include <dune/common/exceptions.hh>

#include <dune/fem/common/coordinate.hh>
#include <dune/fem/function/localfunction/average.hh>
#include <dune/fem/quadrature/geometric/gausslegendre.hh>

#include "continuous.hh"

namespace Dune
{

  namespace Fem
  {

    // Trace
    // -----

    template< class LocalFunction, class Intersection >
    class Trace
    : public ContinuousLocalFunction< typename LocalFunction::FunctionSpaceType, Intersection, Trace< LocalFunction, Intersection > >
    {
      using BaseType = ContinuousLocalFunction< typename LocalFunction::FunctionSpaceType, Intersection, Trace< LocalFunction, Intersection > >;

    public:
      /** \brief local function type */
      using LocalFunctionType = LocalFunction;
      /** \brief intersection type */
      using IntersectionType = Intersection;

      /** \copydoc Dune::Fem::ContinuousLocalFunction::RangeType */
      using RangeType = typename BaseType::RangeType;
      /** \copydoc Dune::Fem::ContinuousLocalFunction::JacobianRangeType */
      using JacobianRangeType = typename BaseType::JacobianRangeType;
      /** \copydoc Dune::Fem::ContinuousLocalFunction::HessianRangeType */
      using HessianRangeType = typename BaseType::HessianRangeType;

    private:
      using LocalGeometryType = typename IntersectionType::LocalGeometry;

    public:
      /** \name Construction
       *  \{
       */
      Trace ( const LocalFunction &localFunction, const IntersectionType &intersection )
        : localFunction_( localFunction ),
          intersection_( intersection ),
          localGeometry_( intersection.geometryInInside() )
      {}

      /** \} */

      /** \name Query methods
       *  \{
       */

      /** \copydoc Dune::Fem::ContinuousLocalFunction::order */
      int order () const { return localFunction().order(); }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::evaluate */
      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        const typename LocalGeometryType::LocalCoordinate y = Dune::Fem::coordinate( x );
        const auto z = localGeometry().global( y );
        localFunction().evaluate( z, value );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::jacobian */
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        DUNE_THROW( NotImplemented, "Method jacobian() not implemented yet" );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::hessian */
      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        DUNE_THROW( NotImplemented, "Method hessian() not implemented yet" );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::evaluateQuadrature */
      template< class Quadrature, class RangeArray >
      void evaluateQuadrature ( const Quadrature &quadrature, RangeArray &values ) const
      {
        const std::size_t nop = quadrature.nop();
        for( std::size_t qp = 0; qp < nop; ++qp )
          evaluate( quadrature[ qp ], values[ qp ] );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::entity */
      const typename BaseType::EntityType &entity () const { return intersection(); }

      /** \} */

      /** \brief Non-interface methods
       *  \{
       */

      /** \brief return local function */
      const LocalFunctionType &localFunction () const { return localFunction_.get(); }

      /** \brief return intersection */
      const IntersectionType &intersection () const { return intersection_.get(); }

      /** \} */

    private:
      const LocalGeometryType &localGeometry () const { return localGeometry_; }

      std::reference_wrapper< const LocalFunction > localFunction_;
      std::reference_wrapper< const Intersection > intersection_;
      LocalGeometryType localGeometry_;
    };



    // LocalAverage
    // ------------

    template< class LocalFunction, class Intersection, class GridPart >
    class LocalAverage< Trace< LocalFunction, Intersection >, GridPart >
    {
    public:
      using LocalFunctionType = Trace< LocalFunction, Intersection >;

      static void apply ( const LocalFunctionType &localFunction, typename LocalFunctionType::RangeType &average )
      {
        using GeometryType = typename LocalFunctionType::EntityType::Geometry;
        const GeometryType &geometry = localFunction.entity().geometry();

        using QuadratureType
          = Dune::Fem::GaussLegendreQuadrature< typename GeometryType::ctype, GeometryType::mydimension >;
        QuadratureType quadrature( geometry.type(), localFunction.order() );

        Dune::Fem::LocalAverageHelper::applyQuadrature( localFunction, geometry, quadrature, average );
      }

      void operator() ( const LocalFunctionType &localFunction, typename LocalFunctionType::RangeType &average ) const
      {
        apply( localFunction, average );
      }
    };



    // trace
    // -----

    template< class LocalFunction, class Intersection >
    Trace< LocalFunction, Intersection >
    trace ( const LocalFunction &localFunction, const Intersection &intersection )
    {
      return Trace< LocalFunction, Intersection >( localFunction, intersection );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_TRACE_HH
