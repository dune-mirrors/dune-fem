#ifndef DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
#define DUNE_FEMPY_PY_DISCRETEFUNCTION_HH

#include <cstddef>

#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

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

      template< class DF >
      inline static std::enable_if_t< std::is_same< decltype( std::declval< const Dune::Fem::DefaultLocalRestrictProlong< typename DF::DiscreteFunctionSpaceType > & >().needCommunication() ), bool >::value >
      registerRestrictProlong ( pybind11::module module, PriorityTag< 1 > )
      {
        typedef typename DF::GridPartType::GridType Grid;

        detail::clsVirtualizedRestrictProlong< Grid >( module ).def( "__init__", [] ( VirtualizedRestrictProlong< Grid > &self, DF &df ) {
            new (&self) VirtualizedRestrictProlong< Grid >( df );
          }, pybind11::keep_alive< 1, 2 >() );
        pybind11::implicitly_convertible< DF, VirtualizedRestrictProlong< Grid > >();
      }

      template< class DF >
      inline static void registerRestrictProlong ( pybind11::module module, PriorityTag< 0 > )
      {}

      template< class DF >
      inline static void registerRestrictProlong ( pybind11::module module )
      {
        registerRestrictProlong< DF >( module, PriorityTag< 42 >() );
      }



      // registerDiscreteFunctionConstructor
      // -----------------------------------

      // specialization for NumPy discrete function, since they require a constructor taking a DoF vector
      template< class Space, class Field, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< Dune::Fem::VectorDiscreteFunction< Space, Dune::FemPy::NumPyVector< Field > >, options... > cls, PriorityTag< 2 > )
      {
        typedef Dune::Fem::VectorDiscreteFunction< Space, Dune::FemPy::NumPyVector< Field > > DF;
        typedef typename DF::VectorType VectorType;

        using pybind11::operator""_a;

        cls.def( pybind11::init( [] ( const Space &space, std::string name, pybind11::buffer dof ) {
            VectorType *vec = new VectorType( std::move( dof ) );
            return new DF( std::move( name ), space, *vec );
          } ), "space"_a, "name"_a, "dof"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 4 >() );
      }

      template< class DF, class... options >
      inline static auto registerDiscreteFunctionConstructor ( pybind11::class_< DF, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_constructible< DF, const std::string &, const typename DF::DiscreteFunctionSpaceType & >::value >
      {
        using pybind11::operator""_a;

        typedef typename DF::DiscreteFunctionSpaceType Space;
        cls.def( "__init__", [] ( DF &self, const Space &space, std::string name ) {
            new (&self) DF( std::move( name ), space );
          }, "space"_a, "name"_a, pybind11::keep_alive< 1, 2 >() );
      }

      template< class DF, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< DF, options... > cls, PriorityTag< 0 > )
      {}

      template< class DF, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< DF, options... > cls )
      {
        registerDiscreteFunctionConstructor( cls, PriorityTag< 42 >() );
      }



      // registerDofVectorBuffer
      // -----------------------

      template < class DofVector, class... options >
      inline static auto registerDofVectorBuffer ( pybind11::class_< DofVector, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_convertible< decltype( std::declval< DofVector >().array().data()[ 0 ] ), typename DofVector::FieldType >::value >
      {
        typedef typename DofVector::FieldType Field;

        cls.def_buffer( [] ( DofVector &self ) -> pybind11::buffer_info {
            return pybind11::buffer_info(
                self.array().data(),                                    /* Pointer to buffer */
                sizeof( Field ),                                        /* Size of one scalar */
                pybind11::format_descriptor< Field >::format(),         /* Python struct-style format descriptor */
                1,                                                      /* Number of dimensions */
                { self.array().size() },                                /* Buffer dimensions */
                { sizeof( Field ) }                                     /* Strides (in bytes) for each index */
            );
          } ); // , pybind11::keep_alive< 0, 1 >() );

        cls.def( "__getitem__", [] ( const DofVector &self, std::size_t index ) -> Field {
            if( index < self.array().size() )
              return self.array().data()[index];
            else
              throw pybind11::index_error();
          });
        cls.def( "__setitem__", [] ( DofVector &self, std::size_t index, Field value ) {
            if( index < self.array().size() )
              return self.array().data()[index] = value;
            else
              throw pybind11::index_error();
          });
      }

      template< class DofVector, class... options >
      inline static void registerDofVectorBuffer ( pybind11::class_< DofVector, options... > cls, PriorityTag< 0 > )
      {}

      template< class DofVector, class... options >
      inline static void registerDofVectorBuffer ( pybind11::class_< DofVector, options... > cls )
      {
        registerDofVectorBuffer( cls, PriorityTag< 42 >() );
      }



      // registerDiscreteFunction
      // ------------------------

      template< class DF, class... options >
      inline static void registerDiscreteFunction ( pybind11::module module, pybind11::class_< DF, options... > cls )
      {
        typedef typename DF::DiscreteFunctionSpaceType Space;
        typedef typename DF::GridPartType GridPart;
        typedef typename DF::RangeType Value;

        detail::registerGridFunction< DF >( module, cls );

        detail::clsVirtualizedGridFunction< GridPart, Value >( module ).def( "__init__", [] ( VirtualizedGridFunction< GridPart, Value > &self, DF &df ) {
            new (&self) VirtualizedGridFunction< GridPart, Value >( pyGridFunction( df ) );
          } );
        pybind11::implicitly_convertible< DF, VirtualizedGridFunction< GridPart, Value > >();

        registerRestrictProlong< DF >( module );

        cls.def_property_readonly( "space", [] ( pybind11::object self ) { return getSpace( self.cast< const DF & >(), self ); } );
        cls.def_property_readonly( "size", [] ( DF &self ) { return self.size(); } );

        registerDiscreteFunctionConstructor( cls );

        cls.def( "copy", [] ( DF &self ) {
            pybind11::object copy = pybind11::cast( new DF( self ), pybind11::return_value_policy::take_ownership );
            // keep alive discrete function space until copy dies, too
            pybind11::detail::keep_alive_impl( copy, pybind11::detail::get_object_handle( &self.space(), pybind11::detail::get_type_info( typeid( Space ) ) ) );
            return copy;
          } );

        cls.def( "clear", [] ( DF &self ) { self.clear(); } );
        cls.def( "assign", [] ( DF &self, const DF &other ) { self.assign( other ); } );

        typedef VirtualizedGridFunction< GridPart, typename Space::RangeType > GridFunction;
        cls.def( "_interpolate", [] ( DF &self, const GridFunction &gf ) {
            Fem::interpolate( gf, self );
          } );
        cls.def( "_interpolate", [] ( DF &self, typename Space::RangeType value ) {
            const auto gf = simpleGridFunction( self.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; }, 0 );
            Fem::interpolate( gf, self );
          } );

        typedef typename DF::DofVectorType DofVector;
        if( !pybind11::already_registered< DofVector >() )
        {
          auto clsDof = pybind11::class_< DofVector >( module, "DofVector", pybind11::buffer_protocol() );

          clsDof.def_property_readonly( "size", [] ( DofVector &self ) { return self.array().size(); } );
          clsDof.def( "__len__", [] ( const DofVector &self ) { return self.array().size(); } );
          clsDof.def( "assign", [] ( DofVector &self, const DofVector &other ) { self = other; } );

          registerDofVectorBuffer( clsDof );
        }

        cls.def_property_readonly( "dofVector", [] ( DF &self ) -> DofVector & { return self.dofVector(); } ); // , pybind11::return_value_policy::reference_internal );
      }

    } // namespace detail



    // registerDiscreteFunction
    // ------------------------

    template< class DF, class... options >
    inline static void registerDiscreteFunction ( pybind11::module module, pybind11::class_< DF, options... > cls )
    {
      detail::registerDiscreteFunction( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
