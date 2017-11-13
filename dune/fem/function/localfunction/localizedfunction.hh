#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALIZEDFUNCTION_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALIZEDFUNCTION_HH

#include <functional>
#include <limits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/common/coordinate.hh>

#include "continuous.hh"

namespace Dune
{

  namespace Fem
  {

    // LocalizedFunction
    // -----------------

    template< class Function, class Entity >
    class LocalizedFunction
    : public DefaultContinuousLocalFunction< Entity, typename Function::FunctionSpaceType::RangeType, LocalizedFunction< Function, Entity > >
    {
      using ThisType = LocalizedFunction< Function, Entity >;
      using BaseType = DefaultContinuousLocalFunction< Entity, typename Function::FunctionSpaceType::RangeType, LocalizedFunction< Function, Entity > >;

    public:
      /** \brief function type */
      using FunctionType = Function;

      /** \copydoc Dune::Fem::GridFunction::DomainType */
      using DomainType = typename BaseType::DomainType;
      /** \copydoc Dune::Fem::GridFunction::RangeType */
      using RangeType = typename BaseType::RangeType;
      /** \copydoc Dune::Fem::GridFunction::JacobianRangeType */
      using JacobianRangeType = typename BaseType::JacobianRangeType;
      /** \copydoc Dune::Fem::GridFunction::HessianRangeType */
      using HessianRangeType = typename BaseType::HessianRangeType;

      /** \copydoc Dune::Fem::GridFunction::EntityType */
      using EntityType = typename BaseType::EntityType;

    private:
      using GeometryType = typename EntityType::Geometry;

    public:
      /** \name Construction
       *  \{
       */

      explicit LocalizedFunction ( const FunctionType &function,
                                   const EntityType &entity,
                                   int order = std::numeric_limits< int >::max() )
        : function_( function ),
          entity_( entity ),
          geometry_( entity.geometry() ),
          order_( order )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      LocalizedFunction ( const ThisType & ) = default;

      /** \brief move constructor */
      LocalizedFunction ( ThisType && ) = default;

      /** \brief assignment operator */
      LocalizedFunction &operator= ( const ThisType & ) = default;

      /** \brief move assignment operator */
      LocalizedFunction &operator= ( ThisType && ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::GridFunction::order */
      int order () const { return order_; }

      /** \copydoc Dune::Fem::GridFunction::evaluate */
      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        const DomainType y = geometry().global( coordinate( x ) );
        function_.evaluate( y, value );
      }

      /** \copydoc Dune::Fem::GridFunction::jacobian */
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        const auto geometry = ThisType::geometry();
        const auto &local = coordinate( x );

        const DomainType global = geometry.global( local );
        function_.jacobian( global, jacobian );

        if( int( EntityType::dimension ) == int( EntityType::Geometry::coorddimension ) )
          return;

        const auto &jacobianTransposed = geometry.jacobianTransposed( local );
        const auto &jacobianInverseTransposed = geometry.jacobianInverseTransposed( local );

        for( int i = 0; i < RangeType::dimension; ++i )
        {
          Dune::FieldVector< typename RangeType::value_type, EntityType::dimension > y;
          jacobianInverseTransposed.mtv( jacobian[ i ], y );
          jacobianTransposed.mtv( y, jacobian[ i ] );
        }
      }

      /** \copydoc Dune::Fem::GridFunction::hessian */
      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        DUNE_THROW( NotImplemented, "Method hessian() not implemented yet" );
      }

      /** \copydoc Dune::Fem::GridFunction::entity */
      const EntityType &entity () const { return entity_.get(); }

    private:
      const GeometryType &geometry () const { return geometry_; }

      FunctionType function_;
      std::reference_wrapper< const EntityType > entity_;
      GeometryType geometry_;
      int order_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALIZEDFUNCTION_HH
