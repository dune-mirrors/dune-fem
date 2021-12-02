#ifndef DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
#define DUNE_FEMPY_PY_DISCRETEFUNCTION_HH


#include <dune/fempy/pybind11/pybind11.hh>

#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>


#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/python/istl/bvector.hh>
#endif // #if HAVE_DUNE_ISTL

#include <cstddef>

#include <string>
#include <type_traits>
#include <utility>
#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#include <dune/fem/io/streams/standardstreams.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fem/common/localcontribution.hh>
// #include <dune/fempy/py/common/numpyvector.hh>
// #include <dune/fempy/py/function/grid.hh>
#include <dune/python/common/numpyvector.hh>
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

      // registerDofVectorBuffer
      // -----------------------
      //register method if data method already available
      template < class DofVector, class... options >
      inline static auto registerDofVectorBuffer ( pybind11::class_< DofVector, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_convertible< decltype( std::declval< DofVector >().array().data()[0] ), typename DofVector::FieldType  >::value >
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
          }); //, pybind11::keep_alive< 0, 1 >() );


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

        cls.def( "__len__", [] ( const DofVector &self ) { return self.array().size(); } );
      }

      template< class DofVector, class... options >
      inline static void registerDofVectorBuffer ( pybind11::class_< DofVector, options... > cls, PriorityTag< 0 > )
      {}

      template< class DofVector, class... options >
      inline static void registerDofVectorBuffer ( pybind11::class_< DofVector, options... > cls )
      {
        registerDofVectorBuffer( cls, PriorityTag< 42 >() );
      }



#if HAVE_DUNE_ISTL
      template< class A , class B>
      inline static const BlockVector<A , B> &getBlockVector (const  BlockVector< A,B > &vector ) noexcept
      {
        return vector;
      }

      template< class A, class B>
      inline static BlockVector<A , B> &setBlockVector (  BlockVector< A,B > &vector ) noexcept
      {
        return vector;
      }

#endif //#if HAVE_DUNE_ISTL

#ifdef PETSC4PY_H // will be set it petsc4py.h was included (so import_petsc4py exists and the python module as well)
      template< class DF , class ... options>
      inline auto addDofVectorBackEnd(pybind11::class_<DF,options...> cls, PriorityTag<3> )
      -> void_t< decltype(std::declval<DF&>().petscVec()) >
      {
        cls.def_property_readonly( "_backend", [] ( DF &self )
        {
          if (import_petsc4py() != 0)
          {                           \
            std::cout << "ERROR: could not import petsc4py\n";
            throw std::runtime_error("Error during import of petsc4py");
          }

          Vec vec = self.dofVector().array();
          pybind11::handle petsc_vec(PyPetscVec_New(vec));
          return petsc_vec;
        }, pybind11::keep_alive<0,1>());
      }
#endif
      template< class DF , class ... options>
      inline auto addDofVectorBackEnd(pybind11::class_<DF,options...> cls, PriorityTag<2> )
      -> void_t< decltype(getBlockVector(std::declval<DF&>().dofVector().array())) >
      {
        typedef typename DF::DofVectorType DofVector;
        //check if BlockVector Is already registered if not register it
        typedef std::decay_t< decltype( getBlockVector( std::declval< DofVector& >().array() ) ) > BlockVector;
        //it's here that I need to add it's name to the type registery
#if HAVE_DUNE_ISTL
        if( !pybind11::already_registered< BlockVector >() )
        {
          Python::registerBlockVector< BlockVector >( cls );
        }
#endif
        cls.def_property_readonly( "_backend", [] ( DF &self )
        -> decltype( getBlockVector(std::declval<DF&>().dofVector().array()) )&
        {
          return self.dofVector().array();
        });
      }
      template< class DF , class ... options>
      inline void addDofVectorBackEnd(pybind11::class_<DF,options...> cls, PriorityTag<1> )
      {
      }
      template< class DF , class ... options>
      inline void addDofVector(pybind11::class_<DF,options...> cls )
      {
        using pybind11::operator""_a;
        typedef typename DF::DofVectorType DofVector;
        if( !pybind11::already_registered< DofVector >() )
        {
          auto clsDof = pybind11::class_< DofVector >( cls, "DofVector", pybind11::buffer_protocol() );
          registerDofVectorBuffer( clsDof );

          clsDof.def_property_readonly( "size", [] ( DofVector &self ) { return self.size(); } );
          clsDof.def( "assign", [] ( DofVector &self, const DofVector &other ) { self = other; }, "other"_a );
          clsDof.def( "scalarProduct", [] ( const DofVector &self, const DofVector &other ) { return self*other; }, "other"_a );
        }
        cls.def_property_readonly( "dofVector", [] ( DF &self ) -> DofVector&
        {
          return self.dofVector();
        });
        addDofVectorBackEnd(cls, PriorityTag<42>());
      }
      // registerRestrictProlong
      // -----------------------

      template< class DF >
      inline static std::enable_if_t< std::is_same< decltype( std::declval< const Dune::Fem::DefaultLocalRestrictProlong< typename DF::DiscreteFunctionSpaceType > & >().needCommunication() ), bool >::value >
      registerRestrictProlong ( pybind11::module module, PriorityTag< 1 > )
      {
        typedef typename DF::GridPartType::GridType Grid;

        detail::clsVirtualizedRestrictProlong< Grid >( module ).def( pybind11::init( [] ( DF &df ) {
            return new VirtualizedRestrictProlong< Grid >( df );
          } ), pybind11::keep_alive< 1, 2 >() );
        pybind11::implicitly_convertible< DF, VirtualizedRestrictProlong< Grid > >();
      }

      template< class DF >
      inline static void registerRestrictProlong ( pybind11::module module, PriorityTag< 0 > )
      {}

      template< class DF >
      inline static void registerRestrictProlong ( pybind11::module module )
      {
        try
        {
          registerRestrictProlong< DF >( module, PriorityTag< 42 >() );
        }
        catch ( const std::invalid_argument& )
        {
          std::cerr << "Warning: Restrict and Prolongation disabled - possibly no HierarchicGrid registered for this GridView " << std::endl;
        }
      }



      // registerDiscreteFunctionConstructor
      // -----------------------------------

      // specialization for NumPy discrete function, since they require a constructor taking a DoF vector
      template< class Space, class Field, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< Dune::Fem::VectorDiscreteFunction< Space, Dune::Python::NumPyVector< Field > >, options... > cls, PriorityTag< 2 > )
      {
        typedef Dune::Fem::VectorDiscreteFunction< Space, Dune::Python::NumPyVector< Field > > DF;
        typedef typename DF::VectorType VectorType;

        using pybind11::operator""_a;

        cls.def( pybind11::init( [] ( const Space &space, std::string name, pybind11::buffer dof ) {
            VectorType *vec = new VectorType( std::move( dof ) );
            DF* df = new DF( std::move( name ), space, *vec );
            // create Python guard object, removing the dof storage once the df disappears
            pybind11::cpp_function remove_dofStorage( [ vec ] ( pybind11::handle weakref ) {
                delete vec;
                weakref.dec_ref();
              } );
            pybind11::weakref weakref( df, remove_dofStorage );
            weakref.release();
            return df;
          } ), "space"_a, "name"_a, "dof"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 4 >() );
      }

      template< class DF, class... options >
      inline static auto registerDiscreteFunctionConstructor ( pybind11::class_< DF, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_constructible< DF, const std::string &, const typename DF::DiscreteFunctionSpaceType & >::value >
      {
        using pybind11::operator""_a;

        cls.def( pybind11::init( [] ( const typename DF::DiscreteFunctionSpaceType &space, std::string name )
          {
            return new DF( std::move( name ), space );
          } ), "space"_a, "name"_a, pybind11::keep_alive< 1, 2 >() );
      }

      template< class DF, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< DF, options... > cls, PriorityTag< 0 > )
      {}

      template< class DF, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< DF, options... > cls )
      {
        registerDiscreteFunctionConstructor( cls, PriorityTag< 42 >() );
      }

      // registerSubDiscreteFunction
      // ---------------------------

      template< class GridFunction, class... options >
      inline static auto registerSubDiscreteFunction ( pybind11::class_< GridFunction, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_same< typename GridFunction::template SubDiscreteFunction< 0 >::Type &, decltype( std::declval< GridFunction & >().template subDiscreteFunction< 0 >() ) >::value >
      {
        cls.def_property_readonly( "components", [] ( pybind11::object self ) -> pybind11::tuple {
            GridFunction &gridFunction = pybind11::cast< GridFunction & >( self );
            pybind11::tuple components( GridFunction::Sequence::size() );
            Hybrid::forEach( typename GridFunction::Sequence(), [ self, &gridFunction, &components ] ( auto i ) {
                assert( pybind11::already_registered< typename GridFunction::template SubDiscreteFunction< i.value >::Type > );
                pybind11::object subFunction = pybind11::cast( &gridFunction.template subDiscreteFunction< i.value >(), pybind11::return_value_policy::reference_internal, self );
                if( subFunction )
                  components[ i.value ] = subFunction;
                else
                  throw pybind11::error_already_set();
              } );
            return components;
          } );
      }

      template< class GridFunction, class... options >
      inline static void registerSubDiscreteFunction ( pybind11::class_< GridFunction, options... > cls, PriorityTag< 0 > )
      {}

      template< class GridFunction, class... options >
      inline static void registerSubDiscreteFunction ( pybind11::class_< GridFunction, options... > cls )
      {
        registerSubDiscreteFunction( cls, PriorityTag< 42 >() );
      }

      template< class GridFunction, class DiscreteFunction, class... options >
      inline static auto registerProjection ( pybind11::class_< DiscreteFunction, options... > cls, PriorityTag< 1 > )
      -> std::enable_if_t<std::is_void< decltype(
           Fem::VtxProjectionImpl::project(
             std::declval<const
                 typename GetGridFunction<
                    typename DiscreteFunction::GridPartType,GridFunction>::value& >(),
             std::declval<DiscreteFunction&>()
           )) >::value>
      {
        using pybind11::operator""_a;
        cls.def( "_project", [] ( DiscreteFunction &self, const GridFunction &gridFunction ) {
            Fem::VtxProjectionImpl::project(
              getGridFunction(self.gridPart(),gridFunction,self.order(),PriorityTag<42>()),
              self);
          }, "gridFunction"_a );
      }

      template< class GridFunction, class DiscreteFunction, class... options >
      inline static void registerProjection ( pybind11::class_< DiscreteFunction, options... > cls, PriorityTag< 0 > )
      {}

      template< class GridFunction, class DiscreteFunction, class... options >
      inline static void registerProjection ( pybind11::class_< DiscreteFunction, options... > cls )
      {
        registerProjection<GridFunction>( cls, PriorityTag< 42 >() );
      }


      // registerDiscreteFunction
      // ------------------------

      template< class DF, class... options >
      inline static void registerDiscreteFunction ( pybind11::module module, pybind11::class_< DF, options... > cls )
      {
        typedef typename DF::DiscreteFunctionSpaceType Space;
        typedef typename DF::GridPartType GridPart;
        typedef typename GridPart::template Codim<0>::EntityType Entity;
        typedef typename DF::RangeType Value;
        const static int dimRange = Value::dimension;

        using pybind11::operator""_a;

        detail::registerGridFunction< DF >( module, cls );

        detail::clsVirtualizedGridFunction< GridPart, Value >( module ).def( pybind11::init( [] ( DF &df ) {
            return new VirtualizedGridFunction< GridPart, Value >( df );
          } ) );
        pybind11::implicitly_convertible< DF, VirtualizedGridFunction< GridPart, Value > >();
#if HAVE_DUNE_VTK
        using GridView = typename DF::GridView;
        using VtkGF = Dune::Vtk::Function<GridView>;
        // register the Function class if not already available
        auto vgfClass = Python::insertClass<VtkGF>(module,"VtkFunction",
            Python::GenerateTypeName("Dune::Vtk::Function",MetaType<GridView>()),
            Python::IncludeFiles{"dune/vtk/function.hh"});
        assert( !vgfClass.second );
        vgfClass.first.def( pybind11::init( [] ( DF &df ) {
            return new VtkGF( localFunction(df), df.name() );
          } ) );
        pybind11::implicitly_convertible<DF,VtkGF>();
#endif

        registerRestrictProlong< DF >( module );

        cls.def_property_readonly( "size", [] ( DF &self ) { return self.size(); } );
        cls.def_property_readonly( "dofsValid", [] ( DF &self ) { return self.dofsValid(); } );
        {
          typedef typename DF::DiscreteFunctionSpaceType DFSpaceType;
          cls.def_property_readonly( "_space", [] ( DF &self ) -> const DFSpaceType& { return self.space(); } );
        }

        registerDiscreteFunctionConstructor( cls );

        cls.def( "copy", [] ( DF &self ) {
            DF *df = new DF( self );
            df->name() = "copy_of_"+self.name();
            pybind11::object copy = pybind11::cast( df, pybind11::return_value_policy::take_ownership );
            // keep alive discrete function space until copy dies, too
            pybind11::detail::keep_alive_impl( copy, pybind11::detail::get_object_handle( &self.space(), pybind11::detail::get_type_info( typeid( Space ) ) ) );
            return copy;
          } );
        cls.def( "copy", [] ( DF &self, const std::string &name ) {
            DF *df = new DF( self );
            df->name() = name;
            pybind11::object copy = pybind11::cast( df , pybind11::return_value_policy::take_ownership );

            // keep alive discrete function space until copy dies, too
            pybind11::detail::keep_alive_impl( copy, pybind11::detail::get_object_handle( &self.space(), pybind11::detail::get_type_info( typeid( Space ) ) ) );
            return copy;
          } );

        cls.def( "clear", [] ( DF &self ) { self.clear(); } );
        // assign is added in __init__ of space
        //cls.def( "assign", [] ( DF &self, const DF &other ) { self.assign( other ); }, "other"_a );
        cls.def( "scalarProductDofs", [] ( DF &self, const DF &other ) { return self.scalarProductDofs( other ); }, "other"_a );
        cls.def( "axpy", [] ( DF &self, double a, const DF &y ) { self.axpy( a, y ); }, "a"_a, "y"_a );
        cls.def("add", [] ( DF &self, const DF &other ) { self += other; }, "other"_a );
        cls.def("sub", [] ( DF &self, const DF &other ) { self -= other; }, "other"_a );
        cls.def("mul", [] ( DF &self, double factor )   { self *= factor; }, "factor"_a );

        cls.def( "_interpolate", [] ( DF &self, typename Space::RangeType value ) {
            const auto gf = simpleGridFunction( self.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; }, 0 );
            Fem::interpolate( gf, self );
          }, "value"_a );
        auto interpol = [&cls](auto *gf) {
          typedef decltype(*gf) GridFunction;
          cls.def( "_interpolate", [] ( DF &self, const GridFunction &gridFunction ) {
              Fem::interpolate( getGridFunction(self.gridPart(),gridFunction,self.order(),PriorityTag<42>()),
                                  self );
            }, "gridFunction"_a );
          registerProjection<GridFunction>(cls);
        };
        addVirtualizedFunctions<GridPart,Value>(interpol);

        if constexpr (dimRange == 1)
        {
          {
            typedef Dune::Python::detail::PyGridFunctionEvaluator< typename GridPart::GridViewType, 0, pybind11::function > LocalEvaluator;
            typedef Dune::Python::SimpleGridFunction< typename GridPart::GridViewType, LocalEvaluator > GridFunction;
            cls.def( "_interpolate", [] ( DF &self, const GridFunction &gridFunction ) {
                Fem::interpolate( simpleGridFunction(self.gridPart(),gridFunction.localEvaluator(),self.order()),
                                  self );
              }, "gridFunction"_a );
            registerProjection<GridFunction>(cls);
          }
          {
            typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;
            typedef std::function<double(const Entity&,const LocalCoordinate&)> stdfct;
            typedef Dune::Python::detail::PyGridFunctionEvaluator< typename GridPart::GridViewType, 0, stdfct > LocalEvaluator;
            typedef Dune::Python::SimpleGridFunction< typename GridPart::GridViewType, LocalEvaluator > GridFunction;
            cls.def( "_interpolate", [] ( DF &self, const GridFunction &gridFunction ) {
                Fem::interpolate( simpleGridFunction(self.gridPart(),gridFunction.localEvaluator(),self.order()),
                                  self );
              }, "gridFunction"_a );
            registerProjection<GridFunction>(cls);
          }
        }

        addDofVector(cls);
        registerSubDiscreteFunction( cls );

        typedef Dune::Fem::AddLocalContribution<DF> AddLocalContrib;
        auto clsAddContrib =
          Dune::Python::insertClass<AddLocalContrib>(module, "AddLocalContribution",
          Dune::Python::GenerateTypeName("AddLocalContribution",MetaType<DF>()),
          Dune::Python::IncludeFiles({"dune/fem/common/localcontribution.hh"}));
        if (clsAddContrib.second)
        {
          auto cls = clsAddContrib.first;
          cls.def("bind",[] (AddLocalContrib &self, const Entity &entity)
              { self.bind(entity); },
            pybind11::keep_alive< 1, 2 >() );
          cls.def("unbind",[] (AddLocalContrib &self)
              { self.unbind(); } );
        }
        cls.def( "addLocalContribution", [] ( DF &self ) {
            auto ret = std::make_unique<AddLocalContrib>(self);
            return ret;
          }, pybind11::keep_alive< 0, 1 >(),
             pybind11::return_value_policy::take_ownership );

        typedef Dune::Fem::SetLocalContribution<DF> SetLocalContrib;
        auto clsSetContrib =
          Dune::Python::insertClass<SetLocalContrib>(module, "SetLocalContribution",
          Dune::Python::GenerateTypeName("SetLocalContribution",MetaType<DF>()),
          Dune::Python::IncludeFiles({"dune/fem/common/localcontribution.hh"}));
        if (clsSetContrib.second)
        {
          auto cls = clsSetContrib.first;
          cls.def("bind",[] (SetLocalContrib &self, const Entity &entity)
              { self.bind(entity); },
            pybind11::keep_alive< 1, 2 >() );
          cls.def("unbind",[] (SetLocalContrib &self)
              { self.unbind(); } );
          cls.def("__getitem__",[] (SetLocalContrib &self, const std::size_t idx)
              { return self[idx]; } );
          cls.def("__setitem__",[] (SetLocalContrib &self, const std::size_t idx, double value)
              { self[idx] = value; } );
        }
        cls.def( "setLocalContribution", [] ( DF &self ) {
            auto ret = std::make_unique<SetLocalContrib>(self);
            return ret;
          }, pybind11::keep_alive< 0, 1 >(),
             pybind11::return_value_policy::take_ownership );

        cls.def( "read", [] ( DF &self, pybind11::bytes state ) {
          std::istringstream stream( state );
          Dune::Fem::StandardInStream inStream(stream);
          self.read( inStream );
          } );
        cls.def( "write", [] ( DF &self ) -> pybind11::bytes {
          std::ostringstream stream;
          Dune::Fem::StandardOutStream outStream(stream);
          self.write( outStream );
          return stream.str();
          } );
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
