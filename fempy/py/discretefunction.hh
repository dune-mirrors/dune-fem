#ifndef DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
#define DUNE_FEMPY_PY_DISCRETEFUNCTION_HH

#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fem/space/common/interpolate.hh>

#include <dune/fempy/py/common/numpyvector.hh>
#include <dune/fempy/py/function/grid.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/restrictprolong.hh>
#include <dune/fempy/py/space.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {
      // registerRestrictProlong
      // -----------------------
      template <class T, class = void>
      struct IsComplete : std::false_type
      {};
      template <class T>
      struct IsComplete< T, decltype(void(sizeof(T))) > : std::true_type
      {};
      template <class Grid, class DF>
      inline static void registerRestrictProlong( pybind11::module module, std::true_type)
      {
        detail::clsVirtualizedRestrictProlong< Grid >( module )
          .def( "__init__", [] ( VirtualizedRestrictProlong< Grid > &instance, DF &df ) {
            new (&instance) VirtualizedRestrictProlong< Grid >( df );
          }, pybind11::keep_alive< 1, 2 >() );
        pybind11::implicitly_convertible< DF, VirtualizedRestrictProlong< Grid > >();
      }
      template <class Grid, class DF>
      inline static void registerRestrictProlong( pybind11::module module, std::false_type)
      {}

      // registerDFConstructor
      // ---------------------

      template< class DF, class... options >
      inline static void registerDFConstructor ( pybind11::class_< DF, options... > cls, std::false_type )
      {}

      template< class DF, class... options >
      inline static void registerDFConstructor ( pybind11::class_< DF, options... > cls, std::true_type )
      {
        using pybind11::operator""_a;

        typedef typename DF::DiscreteFunctionSpaceType Space;
        cls.def( "__init__", [] ( DF &self, const Space &space, std::string name ) {
            new (&self) DF( std::move( name ), space );
          }, "space"_a, "name"_a, pybind11::keep_alive< 1, 2 >() );
      }

      template< class DF, class... options >
      inline static void registerDFConstructor ( pybind11::class_< DF, options... > cls )
      {
        registerDFConstructor( cls, std::is_constructible< DF, const std::string &, const typename DF::DiscreteFunctionSpaceType & >() );
      }


      // registerDofVectorBuffer
      // ----------------

      template < class DofVector, class Cls,
                 typename std::enable_if <
                   std::is_convertible < decltype( std::declval<DofVector>().array().data()[0] ),
                                         typename DofVector::FieldType >::value,
                 int >::type tmp=0 >
      inline static auto registerDofVectorBuffer( Cls cls, int )
      -> decltype(std::declval<DofVector>().array().data(),void())
      {
        typedef typename DofVector::FieldType FieldType;
        cls.def_buffer( [](DofVector &instance) -> pybind11::buffer_info {
            return pybind11::buffer_info(
                instance.array().data(),                /* Pointer to buffer */
                sizeof(FieldType),                      /* Size of one scalar */
                pybind11::format_descriptor<FieldType>::format(), /* Python struct-style format descriptor */
                1,                                      /* Number of dimensions */
                { instance.array().size() },            /* Buffer dimensions */
                { sizeof(FieldType) }                   /* Strides (in bytes) for each index */
            );
          }); // ????  pybind11::keep_alive<0,1>() );
        cls.def( "__getitem__", [] ( const DofVector &self, std::size_t index ) -> FieldType {
            if( index < self.array().size() )
              return self.array().data()[index];
            else
              throw pybind11::index_error();
          });
        cls.def( "__setitem__", [] ( DofVector &self, std::size_t index, FieldType value ) {
            if( index < self.array().size() )
              return self.array().data()[index] = value;
            else
              throw pybind11::index_error();
          });
      }

      template< class DF, class Cls >
      inline static void registerDofVectorBuffer ( Cls cls, long )
      {}

      // registerDiscreteFunction
      // ------------------------

      template< class DF, class Cls >
      inline static void registerDiscreteFunction ( pybind11::module module, Cls cls )
      {
        typedef typename DF::DiscreteFunctionSpaceType Space;
        typedef typename DF::GridPartType GridPart;
        typedef typename DF::RangeType Value;
        typedef typename GridPart::GridType Grid;

        detail::registerGridFunction< DF >( module, cls );

        detail::clsVirtualizedGridFunction< GridPart, Value >( module ).def( "__init__", [] ( VirtualizedGridFunction< GridPart, Value > &instance, DF &df ) {
            new (&instance) VirtualizedGridFunction< GridPart, Value >( pyGridFunction( df ) );
          } );
        pybind11::implicitly_convertible< DF, VirtualizedGridFunction< GridPart, Value > >();

        registerRestrictProlong<Grid,DF>(module, IsComplete<Dune::Fem::DefaultLocalRestrictProlong<Space>>() );

        cls.def_property_readonly( "space", [] ( pybind11::object self ) { return getSpace( self.cast< const DF & >(), self ); } );
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
        cls.def( "_interpolate", [] ( DF &df, const GridFunction &gf ) {
            Fem::interpolate( gf, df );
          } );
        cls.def( "_interpolate", [] ( DF &df, typename Space::RangeType value ) {
            const auto gf = simpleGridFunction( df.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; }, 0 );
            Fem::interpolate( gf, df );
          } );

        typedef typename DF::DofVectorType DofVector;
        if( !pybind11::already_registered< DofVector >() )
        {
          auto clsDof = pybind11::class_< DofVector >( module, "DofVector", pybind11::buffer_protocol() );
          clsDof.def_property_readonly( "size", [] ( DofVector &self ) { return self.array().size(); } );
          clsDof.def( "__len__", [] ( const DofVector &self ) { return self.array().size(); } );
          clsDof.def( "assign", [] ( DofVector &instance, const DofVector &other ) { instance = other; } );
          registerDofVectorBuffer< DofVector >( clsDof, 0 );
        }

        cls.def_property_readonly( "dofVector", [] ( DF &instance ) -> DofVector&{ return instance.dofVector(); }
            ); // ,pybind11::return_value_policy::reference_internal );

        // registerDFBuffer< DF >( cls, 0 );
      }

    } // namespace detail



    // registerDiscreteFunction
    // ------------------------

    template< class DF, class... options >
    inline static void registerDiscreteFunction ( pybind11::module module, pybind11::class_< DF, options... > cls )
    {
      detail::registerDiscreteFunction< DF >( module, cls );
    }

    // special registry for numpy df since they require a constructor
    // taking the dof vector
    template< class Space, class Field, class... options >
    inline static void registerDiscreteFunction ( pybind11::module module,
          pybind11::class_< Dune::Fem::VectorDiscreteFunction< Space, Dune::FemPy::NumPyVector< Field > >, options... > cls )
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
