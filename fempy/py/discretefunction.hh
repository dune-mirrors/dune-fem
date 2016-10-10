#ifndef DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
#define DUNE_FEMPY_PY_DISCRETEFUNCTION_HH

#include <dune/fem/space/common/interpolate.hh>

#include <dune/corepy/pybind11/pybind11.h>
#include <dune/corepy/pybind11/extensions.h>

#include <dune/fempy/py/function/grid.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/restrictprolong.hh>

#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fempy/py/common/numpyvector.hh>

//#if HAVE_EIGEN
//#include <dune/corepy/pybind11/eigen.h>
//#endif // #if HAVE_EIGEN
#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // registerDFConstructor
      // ---------------------

      template< class DF, class Holder, class Alias >
      void registerDFConstructor ( pybind11::class_< DF, Holder, Alias > &cls, std::false_type )
      {}

      template< class DF, class Holder, class Alias >
      void registerDFConstructor ( pybind11::class_< DF, Holder, Alias > &cls, std::true_type )
      {
        using pybind11::operator""_a;

        typedef typename DF::DiscreteFunctionSpaceType Space;
        cls.def( "__init__", [] ( DF &self, const Space &space, std::string name ) {
            new (&self) DF( std::move( name ), space );
          }, "space"_a, "name"_a, pybind11::keep_alive< 1, 2 >() );
      }

      template< class DF, class Holder, class Alias >
      void registerDFConstructor ( pybind11::class_< DF, Holder, Alias > &cls )
      {
        registerDFConstructor( cls, std::is_constructible< DF, const std::string &, const typename DF::DiscreteFunctionSpaceType & >() );
      }



      // registerDFBuffer
      // ----------------

      template <class DF, class Cls>
      auto registerDFBuffer( Cls &cls, int )
      -> decltype(std::declval<DF>().dofVector().array().data(),void())
      {
        typedef typename DF::RangeType Value;
        typedef typename Value::field_type FieldType;
        cls.def_buffer( [](DF &instance) -> pybind11::buffer_info {
            return pybind11::buffer_info(
                instance.dofVector().array().data(),    /* Pointer to buffer */
                sizeof(FieldType),                      /* Size of one scalar */
                pybind11::format_descriptor<FieldType>::format(), /* Python struct-style format descriptor */
                1,                                      /* Number of dimensions */
                { instance.size() },                    /* Buffer dimensions */
                { sizeof(FieldType) }                   /* Strides (in bytes) for each index */
            );
          }); // ????  pybind11::keep_alive<0,1>() );
      }

      template< class DF, class Cls >
      void registerDFBuffer ( Cls &cls, long )
      {}



      // registerDiscreteFunction
      // ------------------------

      template< class DF, class Cls >
      void registerDiscreteFunction ( pybind11::module module, Cls &cls )
      {
        typedef typename DF::DiscreteFunctionSpaceType Space;
        typedef typename DF::GridPartType GridPart;
        typedef typename DF::RangeType Value;
        typedef typename GridPart::GridType Grid;
        typedef typename Value::field_type FieldType;

        detail::registerGridFunction< DF >( module, cls );

        detail::clsVirtualizedGridFunction< GridPart, Value >( module ).def( "__init__", [] ( VirtualizedGridFunction< GridPart, Value > &instance, DF &df ) {
            new (&instance) VirtualizedGridFunction< GridPart, Value >( pyGridFunction( df ) );
          } );
        pybind11::implicitly_convertible< DF, VirtualizedGridFunction< GridPart, Value > >();

        detail::clsVirtualizedRestrictProlong< Grid >( module ).def( "__init__", [] ( VirtualizedRestrictProlong< Grid > &instance, DF &df ) {
            new (&instance) VirtualizedRestrictProlong< Grid >( df );
          }, pybind11::keep_alive< 1, 2 >() );
        pybind11::implicitly_convertible< DF, VirtualizedRestrictProlong< Grid > >();

        cls.def_property_readonly( "space", [] ( DF &df ) -> const typename DF::DiscreteFunctionSpaceType & { return df.space(); } );
        cls.def_property_readonly( "size", [] ( DF &df ) { return df.size(); } );
        cls.def( "clear", [] ( DF &instance ) { instance.clear(); } );

        registerDFConstructor( cls );

        cls.def( "copy", [] ( DF &self ) {
            pybind11::object copy = pybind11::cast( new DF( self ), pybind11::return_value_policy::take_ownership );
            // keep alive discrete function space until copy dies, too
            pybind11::detail::keep_alive_impl( copy, pybind11::detail::get_object_handle( &self.space(), pybind11::detail::get_type_info( typeid( Space ) ) ) );
            return copy;
          } );

        cls.def( "assign", [] ( DF &instance, const DF &other ) { instance.assign(other); } );

        typedef VirtualizedGridFunction< GridPart, typename Space::RangeType > GridFunction;
        cls.def( "interpolate", [] ( DF &df, const GridFunction &gf ) {
            Fem::interpolate( gf, df );
          } );
        cls.def( "interpolate", [] ( DF &df, typename Space::RangeType value ) {
            const auto gf = simpleGridFunction( df.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; }, 0 );
            Fem::interpolate( gf, df );
          } );

        typedef typename DF::DofVectorType DofVector;
        if( !pybind11::already_registered< DofVector >() )
        {
          auto clsDof = pybind11::class_< DofVector >( module, "DofVector" );
        }

        cls.def( "dofVector", [] ( DF &instance ) -> DofVector&{ return instance.dofVector(); },
                 pybind11::return_value_policy::reference_internal );
        cls.def( "assign", [] ( DF &instance, const DofVector &other ) { instance.dofVector() = other; } );

        registerDFBuffer< DF >( cls, 0 );
      }

    } // namespace detail



    // registerDiscreteFunction
    // ------------------------

    template< class DF, class Holder, class Alias >
    void registerDiscreteFunction ( pybind11::module module, pybind11::class_< DF, Holder, Alias > &cls )
    {
      detail::registerDiscreteFunction< DF >( module, cls );
    }

    // special registry for numpy df since they require a constructor
    // taking the dof vector
    template< class Space, class Field, class Holder, class Alias >
    void registerDiscreteFunction ( pybind11::module module,
          pybind11::class_< Dune::Fem::VectorDiscreteFunction< Space, Dune::FemPy::NumPyVector< Field > >, Holder, Alias > &cls )
    {
      using pybind11::operator""_a;

      typedef Dune::Fem::VectorDiscreteFunction<Space,Dune::FemPy::NumPyVector<Field>> DF;
      detail::registerDiscreteFunction<DF>(module,cls);
      typedef typename DF::VectorType VectorType;
      cls.def( "__init__", [] ( DF &self, Space &space, std::string name, pybind11::buffer dof ) {
          VectorType *vec = new VectorType( std::move( dof ) );
          new (&self) DF( std::move( name ), space, *vec );
        }, "space"_a, "name"_a, "dof"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 4 >() );
    }

  } // namespace FemPy

} // namespace Dune


#endif // DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
