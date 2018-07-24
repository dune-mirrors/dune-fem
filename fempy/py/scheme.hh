#ifndef DUNE_FEMPY_PY_SCHEME_HH
#define DUNE_FEMPY_PY_SCHEME_HH

#include <dune/fempy/pybind11/pybind11.hh>

#include <dune/common/typeutilities.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/matrix/colcompspmatrix.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/common/numpyvector.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/pybind11/pybind11.hh>

#if 0
     namespace pybind11
     {
       namespace detail
       {
          template <> class type_caster<_p_Mat>
          {
            public:
            PYBIND11_TYPE_CASTER(Mat, _("mat"));
            // Python to C++
            bool load(handle src, bool)
            {
              value = PyPetscMat_Get(src.ptr());
              return true;
            }
            static handle cast(Mat src, pybind11::return_value_policy policy, handle parent)
            {
               return pybind11::handle(PyPetscMat_New(src));
            }
            operator Mat() { return value; }
          };
        }
      }
#endif

namespace Dune
{

  namespace FemPy
  {

    // registerScheme
    // --------------

    namespace detail
    {

      // registerSchemeConstructor
      // -------------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_constructible< Scheme, const typename Scheme::DiscreteFunctionSpaceType &, typename Scheme::ModelType & >::value >
      {
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Scheme::ModelType ModelType;

        using pybind11::operator""_a;

        cls.def( pybind11::init( [] ( Space &space, ModelType &model ) {
            return new Scheme( space, std::ref(model) );
          } ), "space"_a, "model"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
        cls.def( pybind11::init( [] ( Space &space, ModelType &model, const pybind11::dict &parameters ) {
            return new Scheme( space, std::ref(model), pyParameter( parameters, std::make_shared< std::string >() ) );
          } ), "space"_a, "model"_a, "parameters"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeConstructor( cls, PriorityTag< 42 >() );
      }



      // registerSchemeAssemble
      // ----------------------
      // register assemble method if data method is available (and return value is registered)
      //
      template< class GF, class Scheme, class... options, std::enable_if_t<
            std::is_same< std::decay_t< decltype(std::declval< Scheme >().assemble( std::declval< const GF& >() )) >,
                  typename Scheme::JacobianOperatorType>::value, int > _i = 0 >
      inline static auto registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< const typename Scheme::JacobianOperatorType & >() ) >
      {
        using pybind11::operator""_a;
        cls.def( "assemble", [] ( Scheme &self, const GF &ubar )
          -> const typename Scheme::JacobianOperatorType&
          {
            return self.assemble( ubar );
          }, "ubar"_a, pybind11::return_value_policy::reference_internal );

      }
      template< class GF, class Scheme, class... options >
      inline static void registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {
      }

      template< class Scheme, class... options >
      inline static void registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls )
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        registerSchemeAssemble< DiscreteFunction >( cls, PriorityTag< 42 >() );
        registerSchemeAssemble< VirtualizedGridFunction< GridPart, RangeType > >( cls, PriorityTag< 42 >() );
      }

      template< class Scheme, class... options, std::enable_if_t<
        std::is_constructible<typename Scheme::LinearInverseOperatorType,double,double,Dune::Fem::SolverParameter>::value,int > _i=0 >
      inline static void registerInverseLinearOperator ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def("inverseLinearOperator",[] (Scheme &self, double eps, const pybind11::dict &parameters) {
          return std::make_unique<typename Scheme::LinearInverseOperatorType>
            ( eps, eps, Dune::Fem::SolverParameter(pyParameter( parameters, std::make_shared< std::string >() )) );
        }, "eps"_a, "parameters"_a );
        cls.def("inverseLinearOperator",[] (Scheme &self, typename Scheme::JacobianOperatorType &jOp, double eps, const pybind11::dict &parameters) {
          auto invOp = std::make_unique<typename Scheme::LinearInverseOperatorType>
            ( eps, eps, Dune::Fem::SolverParameter(pyParameter( parameters, std::make_shared< std::string >() )) );
          invOp->bind(jOp);
          return invOp;
        }, "jOp"_a, "eps"_a, "parameters"_a, pybind11::keep_alive<0,2>() );
      }
      template< class Scheme, class... options >
      inline static void registerInverseLinearOperator ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {
      }
      template< class Scheme, class... options >
      inline static void registerInverseLinearOperator ( pybind11::class_< Scheme, options... > cls )
      {
        using pybind11::operator""_a;
        cls.def("inverseLinearOperator",[] (Scheme &self, double eps) {
          return std::make_unique<typename Scheme::LinearInverseOperatorType>( eps, eps );
        }, "eps"_a );
        cls.def("inverseLinearOperator",[] (Scheme &self, typename Scheme::JacobianOperatorType &jOp, double eps) {
          auto invOp = std::make_unique<typename Scheme::LinearInverseOperatorType>( eps, eps );
          invOp->bind(jOp);
          return invOp;
        }, "jOp"_a, "eps"_a, pybind11::keep_alive<0,2>() );
        registerInverseLinearOperator( cls, PriorityTag<42>() );
      }

#if 0
      // registerSchemeGeneralCall
      // -------------------------
      template< class Scheme, class... options >
      inline static auto registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Scheme & >()(
                     std::declval< const VirtualizedGridFunction< typename Scheme::GridPartType, typename Scheme::DiscreteFunctionSpaceType::RangeType > & >(),
                     std::declval< typename Scheme::DiscreteFunctionType & >()
                   ) ) >
      {
        typedef typename Scheme::DiscreteFunctionSpaceType::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        cls.def( "__call__", [] ( Scheme &self, const VirtualizedGridFunction< GridPart, RangeType > &arg, DiscreteFunction &dest ) { self( arg, dest ); } );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeGeneralCall( cls, PriorityTag< 42 >() );
      }

      // registerSchemeModel
      // -------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeModel ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Scheme >().model() ) >
      {
        cls.def_property_readonly( "model", &Scheme::model );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeModel ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeModel ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeModel( cls, PriorityTag< 42 >() );
      }
#endif

      // registerOperatorSetQuadrtureOrder
      // ---------------------------------

      template< class Operator, class... options >
      inline static auto registerOperatorQuadratureOrders ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Operator >().setQuadratureOrders(0,0) ) >
      {
        cls.def( "setQuadratureOrders", &Operator::setQuadratureOrders );
      }

      template< class Operator, class... options >
      inline static void registerOperatorQuadratureOrders ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}

      template< class Operator, class... options >
      inline static void registerOperatorQuadratureOrders ( pybind11::class_< Operator, options... > cls )
      {
        registerOperatorQuadratureOrders( cls, PriorityTag< 42 >() );
      }

      // registerScheme
      // --------------

      template< class Scheme, class... options >
      inline static void registerScheme ( pybind11::module module, pybind11::class_< Scheme, options... > cls )
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;

        using pybind11::operator""_a;

        registerSchemeConstructor( cls );

        cls.def( "_solve", [] ( Scheme &self, const DiscreteFunction &rhs, DiscreteFunction &solution ) {
            auto info = self.solve( rhs, solution );
            return std::map<std::string,std::string> {
                {"converged",std::to_string(info.converged)},
                {"iterations",std::to_string(info.nonlinearIterations)},
                {"linear_iterations",std::to_string(info.linearIterations)}
              };
          } );
        cls.def( "_solve", [] ( Scheme &self, DiscreteFunction &solution ) {
            auto info = self.solve( solution );
            return std::map<std::string,std::string> {
                {"converged",std::to_string(info.converged)},
                {"iterations",std::to_string(info.nonlinearIterations)},
                {"linear_iterations",std::to_string(info.linearIterations)}
              };
          } );

        cls.def( "setErrorMeasure", &Scheme::setErrorMeasure,
                 pybind11::keep_alive<1,2>() );

        cls.def_property_readonly( "dimRange", [] ( Scheme & ) -> int { return DiscreteFunction::FunctionSpaceType::dimRange; } );
        cls.def_property_readonly( "space", [] ( pybind11::object self ) { return detail::getSpace( self.cast< const Scheme & >(), self ); } );
        registerOperatorQuadratureOrders ( cls );

        auto clsInvOp = Dune::Python::insertClass< typename Scheme::LinearInverseOperatorType >
              ( cls, "LinearInverseOperator", Dune::Python::GenerateTypeName(cls,"LinearInverseOperatorType"));
        if( clsInvOp.second )
        {
          Dune::FemPy::detail::registerBasicOperator(clsInvOp.first);
          clsInvOp.first.def("bind",[] (typename Scheme::LinearInverseOperatorType &self,
                                  typename Scheme::JacobianOperatorType &jOp) {
              self.bind(jOp);
          });
          clsInvOp.first.def_property_readonly("iterations",[]( typename Scheme::LinearInverseOperatorType &self) {
              return self.iterations();
          });
        }
        registerInverseLinearOperator( cls );
        registerSchemeAssemble( cls );
        Dune::FemPy::registerOperator(module,cls);
      }

    } // namespace detail

    template< class Scheme, class... options >
    inline static void registerScheme ( pybind11::module module, pybind11::class_< Scheme, options... > cls )
    {
      detail::registerScheme( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SCHEME_HH
