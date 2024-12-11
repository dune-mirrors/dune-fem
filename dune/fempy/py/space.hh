#ifndef DUNE_FEMPY_PY_SPACE_HH
#define DUNE_FEMPY_PY_SPACE_HH

#include <dune/common/hybridutilities.hh>

#include <dune/fem/space/basisfunctionset/codegen.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include <dune/python/common/dynmatrix.hh>
#include <dune/python/common/dynvector.hh>
#include <dune/python/common/fmatrix.hh>
#include <dune/python/common/fvector.hh>

#include <dune/fempy/py/grid/gridpart.hh>

#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // registerSpaceConstructor
      // -------------------------
      template< class Space, class... options >
      void registerSpaceConstructor ( pybind11::class_< Space, options... > cls, std::false_type )
      {}
      template< class Space, class... options >
      void registerSpaceConstructor ( pybind11::class_< Space, options... > cls, std::true_type )
      {
        using pybind11::operator""_a;

        typedef typename Space::GridPartType GridPart;
        typedef typename GridPart::GridViewType GridView;

        cls.def( pybind11::init( [] ( GridView &gridView ) {
            return new Space( gridPart< GridView >( gridView ) );
          }), pybind11::keep_alive< 1, 2 >(), "gridView"_a );
        cls.def(pybind11::pickle(
            [](const pybind11::object &self) { // __getstate__
            Space& spc = self.cast<Space&>();
            pybind11::dict d;
            if (pybind11::hasattr(self, "__dict__")) {
              d = self.attr("__dict__");
              d[pybind11::str("_zero")] = pybind11::none();
            }
            return pybind11::make_tuple(spc.gridPart().gridView(),d);
          },
        [](pybind11::tuple t) { // __setstate__
            if (t.size() != 2)
                throw std::runtime_error("Invalid state in Space::setstate with "+std::to_string(t.size())+" arguments!");
            pybind11::handle gvPtr = t[0];
            /* Create a new C++ instance */
            Space *spc = new Space( gridPart< GridView >(gvPtr) );
            auto py_state = t[1].cast<pybind11::dict>();
            return std::make_pair(spc, py_state);
          }
        ), pybind11::keep_alive< 1, 2 >());
      }

      template< class Space, class... options >
      void registerSpaceConstructorWithOrder ( pybind11::class_< Space, options... > cls, std::false_type )
      {}
      template< class Space, class... options >
      void registerSpaceConstructorWithOrder ( pybind11::class_< Space, options... > cls, std::true_type )
      {
        using pybind11::operator""_a;

        typedef typename Space::GridPartType GridPart;
        typedef typename GridPart::GridViewType GridView;

        cls.def( pybind11::init( [] ( GridView &gridView, int order ) {
            return new Space( gridPart< GridView >( gridView ), order );
          } ), pybind11::keep_alive< 1, 2 >(), "gridView"_a, "order"_a );
        cls.def(pybind11::pickle(
            [](const pybind11::object &self) { // __getstate__
            Space& spc = self.cast<Space&>();
            pybind11::dict d;
            if (pybind11::hasattr(self, "__dict__")) {
              d = self.attr("__dict__");
              d[pybind11::str("_zero")] = pybind11::none();
            }
            return pybind11::make_tuple(spc.gridPart().gridView(), spc.order(), d);
          },
        [](pybind11::tuple t) { // __setstate__
            if (t.size() != 3)
                throw std::runtime_error("Invalid state in Space::setstate with "+std::to_string(t.size())+" arguments!");
            auto order = t[1].cast<int>();
            pybind11::handle gvPtr = t[0];
            /* Create a new C++ instance */
            Space *spc = new Space( gridPart< GridView >(gvPtr), order );
            auto py_state = t[2].cast<pybind11::dict>();
            return std::make_pair(spc, py_state);
          }
        ), pybind11::keep_alive< 1, 2 >());
      }

      template< class Space, class... options >
      void registerSpaceConstructor ( pybind11::class_< Space, options... > cls )
      {
        typedef typename Space::GridPartType GridPart;
        registerSpaceConstructorWithOrder( cls, std::is_constructible< Space, GridPart&, int >() );
        registerSpaceConstructor( cls, std::is_constructible< Space, GridPart& >() );
      }

      template< class Space, class... options >
      inline static auto registerSubSpace ( pybind11::class_< Space, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_same< const typename Space::template SubDiscreteFunctionSpace< 0 >::Type&, decltype( std::declval< Space & >().template subDiscreteFunctionSpace< 0 >() ) >::value >
      {
        cls.def_property_readonly( "components", [] ( pybind11::object self ) -> pybind11::tuple {
            Space &space = pybind11::cast< Space & >( self );
            pybind11::tuple components( Space::Sequence::size() );
            Hybrid::forEach( typename Space::Sequence(), [ self, &space, &components ] ( auto i ) {
                assert( pybind11::already_registered< typename Space::template SubDiscreteFunctionSpace< i.value >::Type > );
                pybind11::object subSpace = pybind11::cast( &space.template subDiscreteFunctionSpace< i.value >(), pybind11::return_value_policy::reference_internal, self );
                if( subSpace )
                  components[ i.value ] = subSpace;
                else
                  throw pybind11::error_already_set();
              } );
            return components;
          } );
      }
      template< class Space, class... options >
      inline static void registerSubSpace ( pybind11::class_< Space, options... > cls, PriorityTag< 0 > )
      {}

      template< class Space, class... options >
      inline static void registerSubSpace ( pybind11::class_< Space, options... > cls )
      {
        registerSubSpace( cls, PriorityTag< 42 >() );
      }

      // registerSpace
      // -------------

      template< class Space, class... options >
      void registerFunctionSpace ( pybind11::handle module, pybind11::class_< Space, options... > cls )
      {
        typedef typename Space::GridPartType GridPart;
        typedef typename GridPart::GridViewType GridView;
        cls.def_property_readonly( "dimRange", [] ( Space & ) -> int { return Space::dimRange; } );
        cls.def_property_readonly( "dimDomain", [] ( Space & ) -> int { return Space::FunctionSpaceType::dimDomain; } );
        cls.def_property_readonly( "grid", [] ( Space &self ) -> const GridView&
        {
          PyErr_WarnEx(PyExc_DeprecationWarning, "attribute 'grid' is deprecated, use 'gridView' instead.", 2);
          return self.gridPart(); }
        );
        cls.def_property_readonly( "gridView", [] ( Space &self ) -> const GridView& { return self.gridPart(); } );
        cls.def_property_readonly( "order", [] ( Space &self ) -> int { return self.order(); } );
        cls.def( "as_ufl", [] ( pybind11::object &self ) -> auto {
              return Dune::FemPy::getSpaceWrapper()(self);
            });
      }

      template< class Space, class... options >
      void registerSpace ( pybind11::handle module, pybind11::class_< Space, options... > cls )
      {
        typedef typename Space::EntityType  EntityType;
        typedef typename Dune::FieldVector< typename Space::DomainFieldType,
                                            Space::GridPartType::dimension >  LocalDomainType;
        typedef typename Space::RangeType   RangeType;
        typedef typename Space::JacobianRangeType  JacobianRangeType;
        typedef typename Space::BlockMapperType    BlockMapperType;
        // this only works when GlobalKeyType is size_t or similar, which it is
        // for all implementations so far
        typedef typename BlockMapperType::GlobalKeyType  GlobalKeyType;

        registerFunctionSpace(module,cls);
        cls.def_property_readonly( "sequence", [] ( Space &self ) -> int { return self.sequence(); } );
        cls.def_property_readonly( "size", [] ( Space &self ) -> int { return self.size(); } );
        cls.def_property_readonly( "primarySize", [] ( Space &self ) -> int { return self.primarySize(); } );
        cls.def_property_readonly( "auxiliarySize", [] ( Space &self ) -> int { return self.auxiliarySize(); } );
        cls.def_property_readonly( "_sizeOfField", [] ( Space &self ) -> unsigned int { return sizeof(typename Space::RangeFieldType); } );
        cls.def_property_readonly( "localBlockSize", [] ( Space &self ) -> unsigned int { return self.localBlockSize; } );

        // operator __len__
        cls.def("__len__", [] ( Space &self ) -> int { return self.size(); } );

        cls.def("localOrder", [] ( Space &self, const EntityType &e) -> int { return self.order(e); } );
        // we can deprecate the following method it's can be replaced by
        // space.mapper(element) which is more flexible
        cls.def("map", [] ( Space &self, const EntityType &e) -> std::vector< GlobalKeyType >
            { std::vector< GlobalKeyType > indices( self.blockMapper().numDofs( e ) );
              // fill vector with dof indices
              self.blockMapper().mapEach(e, [&indices]( const int local, GlobalKeyType global ) { indices[ local ] = global; });
              return indices;
            } );
        cls.def("evaluateBasis", [] ( Space &self, const EntityType &e, const LocalDomainType& xLocal) -> std::vector<RangeType>
            {
              // get basis function set
              const auto basisSet = self.basisFunctionSet( e );
              // evaluate all basis functions at given point xLocal
              std::vector< RangeType > phis( basisSet.size() );
              basisSet.evaluateAll( xLocal, phis );
              return phis;
            } );
        cls.def("jacobianBasis", [] ( Space &spc, const EntityType &e, const LocalDomainType& xLocal) -> std::vector< JacobianRangeType >
            {
              // get basis function set
              const auto basisSet = spc.basisFunctionSet( e );

              std::vector< JacobianRangeType > dPhis( basisSet.size() );
              basisSet.jacobianAll( xLocal, dPhis );
              return dPhis;
            } );
        #if 0 // the return value here is not convertible to Python so needs some export or we transform it directly to numpy
        typedef typename Space::HessianRangeType   HessianRangeType;
        cls.def("hessianBasis", [] ( Space &spc, const EntityType &e, const LocalDomainType& xLocal)
            -> auto // std::vector< HessianRangeType >
            {
              // get basis function set
              const auto basisSet = spc.basisFunctionSet( e );

              std::vector< HessianRangeType > d2Phis_( basisSet.size() );
              basisSet.hessianAll( xLocal, d2Phis_ );

              constexpr int dimR = Space::dimRange;
              constexpr int dimD = Space::FunctionSpaceType::dimDomain;
              typedef Dune::FieldVector<Dune::FieldMatrix<double, dimD, dimD>, dimR> HesRT;
              std::vector< HesRT > d2Phis( basisSet.size() );
              for (std::size_t i=0;i<basisSet.size();++i)
                d2Phis[i] = d2Phis_[i];
              return d2Phis;
            } );
        #endif

        // Note: we need to construct a `NonBlockMapper`.
        // Reason: For example the `BlockedMapper` for `LagrangeSpaces` with
        // different `dimRange` values is always the same type (and the
        // same object even). We can only register this type once. In the
        // registered `__call__` method we then don't have any access
        // anymore to the local block size (`dimRange`) of the space for
        // which the `mapper` was called. To solve this `space.mapper` on
        // the Python side returns a `NonBlockedMapper`.
        typedef Dune::Fem::NonBlockMapper< BlockMapperType, Space::localBlockSize >
                MapperType;
        auto nbMap = Dune::Python::insertClass< MapperType >
              ( module, "NonBlockMapper", Dune::Python::GenerateTypeName(cls,"BlockMapperType"));
        if( nbMap.second )
        {
          auto clsMap = nbMap.first;
          clsMap.def_property_readonly( "localBlockSize", [] ( MapperType &self )
              -> unsigned int { return MapperType::blockSize; } );
          clsMap.def("block", [] ( MapperType &self, const EntityType &e)
              -> std::vector< GlobalKeyType >
              { const BlockMapperType &bm = self.blockMapper();
                std::vector< GlobalKeyType > indices( bm.numDofs( e ) );
                // fill vector with dof indices
                bm.mapEach(e, [&indices]( const int local, GlobalKeyType global )
                    { indices[ local ] = global; });
                return indices;
              } );
          clsMap.def("__call__", [] ( MapperType &self, const EntityType &e)
              -> std::vector< GlobalKeyType >
              {
                std::vector< GlobalKeyType > indices(self.size());
                // fill vector with dof indices
                self.map( e, indices );
                return indices;
              } );
        }
        cls.def( "mapper", [] ( Space &self ) -> auto
            { return MapperType(self.blockMapper()); },
          pybind11::keep_alive<0,1>()
        );

        registerSpaceConstructor( cls );
        registerSubSpace( cls );

        cls.def( "_generateQuadratureCode", []( Space &self,
                 const std::vector<unsigned int> &interiorOrders,
                 const std::vector<unsigned int> &skeletonOrders,
                 const std::string &path,
                 const std::string &filename)
            {
            #ifndef NDEBUG
              std::cout << "Generate code to " << filename << std::endl;
            #endif
              Dune::Fem::generateCode(self, interiorOrders, skeletonOrders, path, filename);
            } );

      }

    } // namespace detail



    // registerSpace
    // -------------

    template< class Space, class... options >
    void registerSpace ( pybind11::handle module, pybind11::class_< Space, options... > cls )
    {
      detail::registerSpace( module, cls );
    }
  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
