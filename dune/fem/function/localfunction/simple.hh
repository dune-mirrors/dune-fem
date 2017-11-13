#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_SIMPLE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_SIMPLE_HH

#include <functional>
#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/fem/common/coordinate.hh>
#include <dune/fem/function/localfunction/continuous.hh>

namespace Dune
{

  namespace Fem
  {

    // SimpleLocalFunction
    // -------------------

    template< class Entity, class Function >
    class SimpleLocalFunction
      : public DefaultContinuousLocalFunction< Entity, typename std::result_of< Function( Entity, typename Entity::Geometry::LocalCoordinate ) >::type, SimpleLocalFunction< Entity, Function > >
    {
      using BaseType = DefaultContinuousLocalFunction< Entity, typename std::result_of< Function( Entity, typename Entity::Geometry::LocalCoordinate ) >::type, SimpleLocalFunction< Entity, Function > >;

    public:
      /** \copydoc Dune::Fem::ContinuousLocalFunction::EntityType */
      using EntityType = typename BaseType::EntityType;

      /** \name Construction
       *  \{
       */

      SimpleLocalFunction ( const EntityType &entity, Function function, int order )
        : entity_( entity ),
          function_( std::move( function ) ),
          order_( order )
      {}

      /** \} */

      /** \copydoc Dune::Fem::ContinuousLocalFunction::order */
      int order () const
      {
        return order_;
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::evaluate */
      template< class Point >
      void evaluate ( const Point &x, typename BaseType::RangeType &value ) const
      {
        value = function_( entity(), coordinate( x ) );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::jacobian */
      template< class Point >
      void jacobian ( const Point &x, typename BaseType::JacobianRangeType &jacobian ) const
      {
        DUNE_THROW( NotImplemented, "Method jacobian() not implemented" );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::hessian */
      template< class Point >
      void hessian ( const Point &x, typename BaseType::HessianRangeType &hessian ) const
      {
        DUNE_THROW( NotImplemented, "Method hessian() not implemented" );
      }

      /** \copydoc Dune::Fem::ContinuousLocalFunction::entity */
      const EntityType &entity () const
      {
        return entity_.get();
      }

    private:
      std::reference_wrapper< const EntityType > entity_;
      Function function_;
      int order_;
    };



    // simpleLocalFunction
    // -------------------

    template< class Entity, class Function >
    SimpleLocalFunction< Entity, Function >
    simpleLocalFunction ( const Entity &entity, Function function, int order )
    {
      return SimpleLocalFunction< Entity, Function >( entity, std::move( function ), order );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_SIMPLE_HH
