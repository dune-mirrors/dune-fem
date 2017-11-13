#ifndef DUNE_FEM_FUNCTION_COMMON_GRIDFUNCTION_HH
#define DUNE_FEM_FUNCTION_COMMON_GRIDFUNCTION_HH

#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/function.hh>

namespace Dune
{

  namespace Fem
  {

    // GridFunctionInterface
    // ---------------------

    template< class GridPart, class LocalFunction, class Implementation >
    class GridFunction
    : public Dune::Fem::Function< typename LocalFunction::FunctionSpaceType, Implementation >,
      public Dune::Fem::HasLocalFunction
    {
      using BaseType = Dune::Fem::Function< typename LocalFunction::FunctionSpaceType, Implementation >;

    public:
      /** \brief domain type */
      using DomainType = typename BaseType::DomainType;
      /** \brief range type */
      using RangeType = typename BaseType::RangeType;
      /** \brief jacobian range type */
      using JacobianRangeType = typename BaseType::JacobianRangeType;
      /** \brief hessian range type */
      using HessianRangeType = typename BaseType::HessianRangeType;

      /** \brief grid part type */
      using GridPartType = GridPart;
      /** \brief local function type */
      using LocalFunctionType = LocalFunction;

      /** \brief entity type */
      using EntityType = typename LocalFunctionType::EntityType;

    protected:
      GridFunction () = default;

    public:
      /** \name Evaluation in global coordinates
       *  \{
       */

      /** \brief evaluate grid function */
      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().evaluate( x, value ) );
        impl().evaluate( x, value );
      }

      /** \brief evaluate jacobian of grid function */
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().jacobian( x, jacobian ) );
        impl().jacobian( x, jacobian );
      }

      /** \brief evaluate hessian of grid function */
      void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().hessian( x, hessian ) );
        impl().hessian( x, hessian );
      }

      /** \} */

      /** \name Evaluation in local coordinates
       *  \{
       */

      /** \brief return local function */
      LocalFunctionType localFunction ( const EntityType &entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().localFunction( entity ) );
        return impl().localFunction( entity );
      }

      /** \brief return grid part */
      const GridPartType &gridPart () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().gridPart() );
        return impl().gridPart();
      }

      /** \} */

    private:
      const Implementation &impl () const
      {
        return static_cast< const Implementation & >( *this );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_GRIDFUNCTION_HH
