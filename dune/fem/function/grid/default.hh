#ifndef DUNE_FEM_FUNCTION_GRID_DEFAULT_HH
#define DUNE_FEM_FUNCTION_GRID_DEFAULT_HH

#include <functional>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/common/entitysearch.hh>

#include "gridfunction.hh"

namespace Dune
{

  namespace Fem
  {

    // DefaultGridFunction
    // -------------------

    template< class GridPart, class LocalFunction, class Implementation >
    class DefaultGridFunction
    : public Dune::Fem::GridFunction< GridPart, LocalFunction, Implementation >
    {
      using ThisType = DefaultGridFunction< GridPart, LocalFunction, Implementation >;
      using BaseType = Dune::Fem::GridFunction< GridPart, LocalFunction, Implementation >;

    public:
      /** \copydoc Dune::Fem::GridFunction::DomainType */
      using DomainType = typename BaseType::DomainType;
      /** \copydoc Dune::Fem::GridFunction::RangeType */
      using RangeType = typename BaseType::RangeType;
      /** \copydoc Dune::Fem::GridFunction::JacobianRangeType */
      using JacobianRangeType = typename BaseType::JacobianRangeType;
      /** \copydoc Dune::Fem::GridFunction::HessianRangeType */
      using HessianRangeType = typename BaseType::HessianRangeType;

      /** \copydoc Dune::Fem::GridFunction::GridPartType */
      using GridPartType = typename BaseType::GridPartType;
      /** \copydoc Dune::Fem::GridFunction::EntityType */
      using EntityType = typename BaseType::EntityType;

    private:
      using EntitySearchType = Dune::Fem::EntitySearch< GridPartType, EntityType::codimension >;

    public:
      using BaseType::localFunction;

      /** \name Public member methods
       *  \{
       */

      explicit DefaultGridFunction ( const GridPartType &gridPart )
        : gridPart_( gridPart ),
          entitySearch_( std::make_shared< EntitySearchType >( gridPart ) )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      DefaultGridFunction ( const ThisType & ) = default;

      /** \brief move constructor */
      DefaultGridFunction ( ThisType && ) = default;

      /** \brief assignment operator */
      DefaultGridFunction &operator= ( const ThisType & ) = default;

      /** \brief move assignment operator */
      DefaultGridFunction &operator= ( ThisType && ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::GridFunction::evaluate */
      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        EntityType entity = (*entitySearch_)( x );
        const auto localFunction = this->localFunction( entity );
        const auto geometry = entity.geometry();
        localFunction.evaluate( geometry.local( x ), value );
      }

      /** \copydoc Dune::Fem::GridFunction::jacobian */
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        DUNE_THROW( NotImplemented, "Method jacobian() not implemented yet" );
      }

      /** \copydoc Dune::Fem::GridFunction::hessian */
      void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        DUNE_THROW( NotImplemented, "Method hessian() not implemented yet" );
      }

      /** \copydoc Dune::Fem::GridFunction::gridPart */
      const GridPart &gridPart () const { return gridPart_.get(); }

      /** \} */

    private:
      std::reference_wrapper< const GridPartType > gridPart_;
      std::shared_ptr< EntitySearchType > entitySearch_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_GRID_DEFAULT_HH
