#ifndef DUNE_FEMPY_PYGRIDFUNCTION_HH
#define DUNE_FEMPY_PYGRIDFUNCTION_HH

#include <string>
#include <tuple>

#include <dune/common/std/utility.hh>

#include <dune/fempy/gridfunction.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    template< class GridFunction >
    struct RegisterGridFunction;

    template< class GridPart, int dimR >
    struct RegisterGridFunction< GridFunction< GridPart, dimR > >
    {
      typedef GridFunction< GridPart, dimR > GF;
      static void apply ( pybind11::module module, const std::string &name = "GridFunction" )
      {
        static const std::string clsName = name;
        pybind11::class_< GF, std::shared_ptr< GF > > gf( module, clsName.c_str() );
        gf.def( "getDimRange", &GridFunction< GridPart, dimR >::getDimRange );
      }
    };

    template< class GF >
    struct RegisterGridFunction
    {
      typedef GridFunction< typename GF::GridPartType, GF::dimRange > Base;

      static pybind11::class_< GF > apply ( pybind11::module module, const std::string &name )
      {
        static const std::string clsName = name;
        pybind11::class_< GF, std::shared_ptr< GF > > gf( module, clsName.c_str(), pybind11::base< Base >() );
        gf.def( "getDimRange", &GF::getDimRange );
        return gf;
      }

      template< class Init >
      static pybind11::class_< GF > apply ( pybind11::module module, const std::string &name, Init init )
      {
        return apply( module, name ).gf.def( init );
      }
    };

    template< template< int > class GridFunction >
    struct RegisterGridFunctions
    {
      template< int... i, class... Args >
      static void apply ( pybind11::module module, Std::integer_sequence< int, i... >, const std::string &name, Args &&... args )
      {
        std::ignore = std::make_tuple( (RegisterGridFunction< GridFunction< i > >::apply( module, name + std::to_string( i ), std::forward< Args >( args )... ), i)... );
      }
    };

    template< template< class, int, class > class GF, class GridPart, class RF = double >
    struct GridFunctionsBinder
    {
      template< int dimR >
      using GridFunction = GF< GridPart, dimR, RF >;
    };

    template< class GridFunction, class... Args >
    void registerGridFunction ( Args &&... args )
    {
      RegisterGridFunction< GridFunction >::apply( std::forward< Args >( args )... );
    }

    template< class GridPart, int... i >
    void registerGridFunctionInterface ( pybind11::module module, Std::integer_sequence< int, i... >, const std::string &name = "GridFunction" )
    {
      typedef GridFunctionsBinder< GridFunction, GridPart > Binder;
      RegisterGridFunctions< Binder::template GridFunction >::apply( module, Std::integer_sequence< int, i... >(), name );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PYGRIDFUNCTION_HH
