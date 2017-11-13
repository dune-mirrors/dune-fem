#ifndef DUNE_FEM_FUNCTION_GRID_LOCALIZABLEFUNCTION_HH
#define DUNE_FEM_FUNCTION_GRID_LOCALIZABLEFUNCTION_HH

#include <limits>
#include <string>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

#include <dune/fem/function/localfunction/localizedfunction.hh>

#include "gridfunction.hh"

namespace Dune
{

  namespace Fem
  {

    // LocalizableFunction
    // -------------------

    template< class Function, class GridPart >
    class LocalizableFunction
    : public Dune::Fem::GridFunction< GridPart, LocalizedFunction< Function, typename GridPart::template Codim< 0 >::EntityType >, LocalizableFunction< Function, GridPart > >
    {
      using ThisType = LocalizableFunction< Function, GridPart >;
      using BaseType = Dune::Fem::GridFunction< GridPart, LocalizedFunction< Function, typename GridPart::template Codim< 0 >::EntityType >, LocalizableFunction< Function, GridPart > >;

    public:
      /** \brief function type */
      using FunctionType = Function;

      /** \copydoc Dune::Fem::GridFunction::GridPartType */
      using GridPartType = typename BaseType::GridPartType;
      /** \copydoc Dune::Fem::GridFunction::EntityType */
      using EntityType = typename BaseType::EntityType;
      /** \copydoc Dune::Fem::GridFunction::LocalFunctionType */
      using LocalFunctionType = typename BaseType::LocalFunctionType;

      /** \copydoc Dune::Fem::GridFunction::FunctionSpaceType */
      using FunctionSpaceType = typename BaseType::FunctionSpaceType;
      /** \copydoc Dune::Fem::GridFunction::DomainType */
      using DomainType = typename BaseType::DomainType;
      /** \copydoc Dune::Fem::GridFunction::RangeType */
      using RangeType = typename BaseType::RangeType;
      /** \copydoc Dune::Fem::GridFunction::JacobianRangeType */
      using JacobianRangeType = typename BaseType::JacobianRangeType;
      /** \copydoc Dune::Fem::GridFunction::HessianRangeType */
      using HessianRangeType = typename BaseType::HessianRangeType;

      /** \brief discrete function space adapter */
      using DiscreteFunctionSpaceType = Dune::Fem::DiscreteFunctionSpaceAdapter< FunctionSpaceType, GridPartType >;

      /** \name Construction
       *  \{
       */

      LocalizableFunction ( const std::string &name,
                            const FunctionType &function,
                            const GridPartType &gridPart,
                            int order = std::numeric_limits< int >::max() )
        : name_( name ),
          function_( function ),
          space_( gridPart, order )
      {}

      LocalizableFunction ( const std::string &name,
                            const GridPartType &gridPart,
                            int order = std::numeric_limits< int >::max() )
        : LocalizableFunction( name, FunctionType(), gridPart, order )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      LocalizableFunction ( const ThisType & ) = default;

      /** \brief move constructor */
      LocalizableFunction ( ThisType && ) = default;

      /** \brief assignment operator */
      ThisType &operator= ( const ThisType & ) = default;

      /** \brief move assignment operator */
      ThisType &operator= ( ThisType && ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::GridFunction::gridPart */
      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        function().evaluate( x, value );
      }

      /** \copydoc Dune::Fem::GridFunction::gridPart */
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        function().jacobian( x, jacobian );
      }

      /** \copydoc Dune::Fem::GridFunction::gridPart */
      void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        function().hessian( x, hessian );
      }

      /** \copydoc Dune::Fem::GridFunction::gridPart */
      LocalFunctionType localFunction ( const EntityType &entity ) const
      {
        return LocalFunctionType( function_, entity, space().order() );
      }

      /** \copydoc Dune::Fem::GridFunction::gridPart */
      const GridPartType &gridPart () const { return space().gridPart(); }

      /** \} */

      /** \name Non-interface methods
       *  \{
       */

      /** \brief return discrete function space adapter */
      const DiscreteFunctionSpaceType &space () const { return space_; }

      /** \brief return name */
      const std::string &name () const { return name_; }

      /** \brief return function */
      FunctionType &function () { return function_; }

      /** \brief return function */
      const FunctionType &function () const { return function_; }

      /** \} */

    private:
      std::string name_;
      FunctionType function_;
      DiscreteFunctionSpaceType space_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_GRID_LOCALIZABLEFUNCTION_HH
