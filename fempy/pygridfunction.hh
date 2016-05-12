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

    // registerLocalFunction
    // ---------------------

    template< class LocalFunction >
    pybind11::class_< LocalFunction > registerLocalFunction ( pybind11::handle scope, const char *clsName = "LocalFunction" )
    {
      typedef typename LocalFunction::LocalCoordinateType LocalCoordinate;

      pybind11::class_< LocalFunction > cls( scope, clsName );

      cls.def_property_readonly( "dimRange", [] ( LocalFunction & ) -> int { return LocalFunction::RangeType::dimension; } );
      cls.def( "evaluate", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::RangeType value;
          lf.evaluate( x, value );
          return value;
        } );
      cls.def( "jacobian", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::JacobianRangeType jacobian;
          lf.jacobian( x, jacobian );
          return jacobian;
        } );

      return cls;
    }



    // registerGridFunction
    // --------------------

    template< class GridFunction >
    pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
    {
      typedef typename GridFunction::LocalFunctionType LocalFunction;
      typedef typename LocalFunction::EntityType Entity;

      pybind11::class_< GridFunction > cls( scope, clsName );

      registerLocalFunction< LocalFunction >( cls );

      cls.def( "__repr__", [] ( GridFunction &gf ) -> std::string {
          return "GridFunction< " + std::to_string( GridFunction::RangeType::dimension ) + " >(name = " + gf.name() + ")";
        } );

      cls.def_property_readonly( "dimRange", [] ( GridFunction &gf ) -> int { return GridFunction::RangeType::dimension; } );

      cls.def_property_readonly( "name", [] ( GridFunction &gf ) -> std::string { return gf.name(); } );

      cls.def( "localFunction", [] ( const GridFunction &gf, const Entity &entity ) -> LocalFunction {
          return gf.localFunction( entity );
        }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >() );

      return cls;
    }



    // PyFunction
    // ----------

    template< class Domain, class Range >
    struct PyFunction;

    template< class DF, int dimD, class RF, int dimR >
    struct PyFunction< FieldVector< DF, dimD >, FieldVector< RF, dimR > >
    {
      typedef Fem::FunctionSpace< DF, RF, dimD, dimR > FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      explicit PyFunction ( pybind11::function evaluate )
        : evaluate_( std::move( evaluate ) )
      {}

      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        pybind11::gil_scoped_acquire acq;
        pybind11::object v( evaluate_.call( x ) );
        value = v.template cast< RangeType >();
      }

      void jacobian ( const DomainType &x, JacobianRangeType &value ) const
      {
        DUNE_THROW( NotImplemented, "PyFunction::jacobian not implemented, yet." );
      }

    private:
      pybind11::function evaluate_;
    };

    /***********************************************************************/

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

    template< class GridPart, int... i >
    void registerGridFunctionInterface ( pybind11::module module, Std::integer_sequence< int, i... >, const std::string &name = "GridFunction" )
    {
      typedef GridFunctionsBinder< GridFunction, GridPart > Binder;
      RegisterGridFunctions< Binder::template GridFunction >::apply( module, Std::integer_sequence< int, i... >(), name );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PYGRIDFUNCTION_HH
