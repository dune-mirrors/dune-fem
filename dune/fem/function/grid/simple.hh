#ifndef DUNE_FEM_FUNCTION_GRID_SIMPLE_HH
#define DUNE_FEM_FUNCTION_GRID_SIMPLE_HH

#include <string>
#include <utility>

#include <dune/fem/space/common/discretefunctionspace.hh>

#include <dune/fem/function/grid/default.hh>
#include <dune/fem/function/localfunction/simple.hh>

namespace Dune
{

  namespace Fem
  {

    // SimpleGridFunction
    // ------------------

    template< class GridPart, class Function >
    class SimpleGridFunction
      : public DefaultGridFunction< GridPart, SimpleLocalFunction< typename GridPart::template Codim< 0 >::EntityType, Function >, SimpleGridFunction< GridPart, Function > >
    {
      using BaseType = DefaultGridFunction< GridPart, SimpleLocalFunction< typename GridPart::template Codim< 0 >::EntityType, Function >, SimpleGridFunction< GridPart, Function > >;

    public:
      /** \copydoc Dune::Fem::GridFunction::GridPartType */
      using GridPartType = typename BaseType::GridPartType;
      /** \copydoc Dune::Fem::GridFunction::FunctionSpaceType */
      using FunctionSpaceType = typename BaseType::FunctionSpaceType;
      /** \brief discrete function space adapter */
      using DiscreteFunctionSpaceType = Dune::Fem::DiscreteFunctionSpaceAdapter< FunctionSpaceType, GridPartType >;

      /** \copydoc Dune::Fem::GridFunction::EntityType */
      using EntityType = typename BaseType::EntityType;
      /** \copydoc Dune::Fem::GridFunction::LocalFunctionType */
      using LocalFunctionType = typename BaseType::LocalFunctionType;

      /** \name Construction
       *  \{
       */

      SimpleGridFunction ( std::string name, const GridPartType &gridPart, Function function, int order )
        : BaseType( gridPart ),
          name_( std::move( name ) ),
          function_( std::move( function ) ),
          space_( gridPart, order )
      {}

      /** \} */

      /** \copydoc Dune::Fem::GridFunction::localFunction */
      LocalFunctionType localFunction ( const EntityType &entity ) const noexcept
      {
        return simpleLocalFunction( entity, function_, space().order() );
      }

      /** \brief return name */
      const std::string &name () const noexcept
      {
        return name_;
      }

      /** \brief return discrete function space adapter */
      const DiscreteFunctionSpaceType &space () const noexcept
      {
        return space_;
      }

    private:
      std::string name_;
      Function function_;
      DiscreteFunctionSpaceType space_;
    };



    // simpleGridFunction
    // -------------------

    template< class GridPart, class Function >
    SimpleGridFunction< GridPart, Function >
    simpleGridFunction ( std::string name, const GridPart &gridPart, Function function, int order )
    {
      return SimpleGridFunction< GridPart, Function >( std::move( name ), gridPart, std::move( function ), order );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_GRID_SIMPLE_HH
