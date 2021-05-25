#ifndef DUNE_FEMPY_PY_OPERATOR_HH
#define DUNE_FEMPY_PY_OPERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fempy/parameter.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/pybind11/pybind11.hh>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/stl_bind.h>

#if HAVE_DUNE_ISTL
#include <dune/istl/bcrsmatrix.hh>
#include <dune/python/istl/bcrsmatrix.hh>
#endif // #if HAVE_DUNE_ISTL

namespace Dune
{

  namespace FemPy
  {
    template< class Operator, class... options >
    inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls );

    namespace detail
    {

      // GeneralGridFunction
      // -------------------

      template< class DiscreteFunction >
      using GeneralGridFunction = VirtualizedGridFunction< typename DiscreteFunction::GridPartType, typename DiscreteFunction::RangeType >;



      // registerGeneralOperatorCall
      // ---------------------------
      template< class Operator, class... options, decltype( std::declval< const Operator & >()( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::RangeFunctionType & >() ), 0 ) = 0 >
      inline static void registerGeneralOperatorCall ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "__call__", [] ( Operator &self, const GeneralGridFunction< typename Operator::DomainFunctionType > &u, typename Operator::RangeFunctionType &w ) { self( u, w ); }, "u"_a, "w"_a );
      }
      template< class Operator, class... options >
      inline static void registerGeneralOperatorCall ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}



      // registerGeneralOperatorJacobian
      // -------------------------------

      template< class Operator, class... options, decltype( std::declval< const Operator & >().jacobian( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::JacobianOperatorType& >() ), 0 ) = 0 >
      inline static void registerGeneralOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "jacobian", [] ( Operator &self, const GeneralGridFunction< typename Operator::DomainFunctionType > &u, typename Operator::JacobianOperatorType &jOp ) { self.jacobian( u, jOp ); jOp.finalize(); }, "u"_a, "jOp"_a );
      }
      template< class Operator, class... options >
      inline static void registerGeneralOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}

      // registerOperatorJacobian
      // ------------------------

      template< class Operator, class... options, decltype( std::declval< const Operator & >().jacobian( std::declval< const typename Operator::DomainFunctionType & >(), std::declval< typename Operator::JacobianOperatorType& >() ), 0 ) = 0 >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;

        cls.def( "jacobian", [] ( Operator &self, const typename Operator::DomainFunctionType &u, typename Operator::JacobianOperatorType &jOp ) { self.jacobian( u, jOp ); jOp.finalize(); }, "u"_a, "jOp"_a );
      }

      template< class Operator, class... options >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}

      // registerConstraints
      // -------------------
      template < class Op, class DF, typename = void >
      struct AddDirichletBC
      { static constexpr bool value = false; };
      template < class Op, class DF>
      struct AddDirichletBC<Op,DF,std::enable_if_t<std::is_void< decltype( std::declval<const Op>().
                  setConstraints( std::declval<DF&>() ) )>::value > >
      { static constexpr bool value = !std::is_same_v<typename Op::DirichletBlockVector,void>; };
      template< class Operator, class... options >
      inline static void registerOperatorConstraints ( pybind11::class_< Operator, options... > cls )
      {
        typedef typename Operator::DomainFunctionType DomainFunction;
        typedef typename Operator::RangeFunctionType  RangeFunction;
        if constexpr (AddDirichletBC<Operator,DomainFunction>::value)
        {
          cls.def( "setConstraints", [] ( Operator &self, DomainFunction &u) { self.setConstraints( u ); } );
          cls.def( "setConstraints", [] ( Operator &self, const typename DomainFunction::RangeType &value, DomainFunction &u) { self.setConstraints( value, u ); } );
          cls.def( "setConstraints", [] ( Operator &self, const DomainFunction &u, RangeFunction &v) { self.setConstraints( u,v ); } );
          cls.def( "subConstraints", [] ( Operator &self, const DomainFunction &u, RangeFunction &v) { self.subConstraints( u,v ); } );
          using DirichletBlockVector = typename Operator::DirichletBlockVector;
          pybind11::bind_vector<DirichletBlockVector>(cls, "DirichletBlockVector");
          cls.def_property_readonly( "dirichletBlocks",  [] ( Operator &self ) -> auto& { return self.dirichletBlocks(); } );
        }
      }

      template< class Operator, class... options >
      inline static auto registerBasicOperator ( pybind11::class_< Operator, options... > cls )
      {
        typedef typename Operator::DomainFunctionType DomainFunction;
        typedef typename Operator::RangeFunctionType RangeFunction;

        using pybind11::operator""_a;

        cls.def( "__call__", [] ( Operator &self, const DomainFunction &u, RangeFunction &w ) { self( u, w ); }, "u"_a, "w"_a );
        registerGeneralOperatorCall( cls, PriorityTag< 42 >() );
      }

      // registerOperatorSetQuadrtureOrder
      // ---------------------------------

      template< class Operator, class... options >
      inline static auto registerOperatorQuadratureOrders ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Operator >().setQuadratureOrders(0,0) ) >
      {
        cls.def( "setQuadratureOrders", [](Operator &self, unsigned int interior, unsigned int surface)
            { self.setQuadratureOrders(interior,surface); } );
      }

      template< class Operator, class... options >
      inline static void registerOperatorQuadratureOrders ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}

      template< class Operator, class... options >
      inline static void registerOperatorQuadratureOrders ( pybind11::class_< Operator, options... > cls )
      {
        registerOperatorQuadratureOrders( cls, PriorityTag< 42 >() );
      }

      template< class Operator, class... options,
        decltype( std::declval< const Operator & >().domainSpace(), 0 ) = 0 >
      inline static void registerOperatorSpaces ( pybind11::class_< Operator, options... > cls, PriorityTag<1> )
      {
        cls.def_property_readonly( "domainSpace", [] ( Operator &self ) -> auto& { return self.domainSpace(); } );
        cls.def_property_readonly( "rangeSpace", [] ( Operator &self) -> auto& { return self.rangeSpace(); } );
      }
      template< class Operator, class... options >
      inline static void registerOperatorSpaces ( pybind11::class_< Operator, options... > cls, PriorityTag<0> )
      {}

      // registerOperator
      // ----------------

      template< class Operator, class... options >
      inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
      {
        using pybind11::operator""_a;

        registerBasicOperator(cls);

        registerOperatorJacobian( cls, PriorityTag< 42 >() );
        registerGeneralOperatorJacobian( cls, PriorityTag< 42 >() );
        registerOperatorConstraints( cls );
        registerOperatorQuadratureOrders ( cls );

        registerOperatorSpaces( cls, PriorityTag< 42 >() );
      }

      //////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////

      // AddMatrixBackend
      // ----------------

      template< class AssembledLinearOperator, class... options>
      inline static void addMatrixBackend ( pybind11::class_< AssembledLinearOperator, options... > cls, PriorityTag< 0 > )
      {
      }

#ifdef PETSC4PY_H // will be set it petsc4py.h was included (so import_petsc4py exists and the python module as well)
      template< class Mat >
      inline static const Mat &getPetscMatrix ( const Mat &matrix ) noexcept
      {
        return matrix;
      }
      template< class Operator, class... options>
      inline static auto addMatrixBackend ( pybind11::class_< Operator, options... > cls, PriorityTag< 4 > )
      -> void_t<decltype( getPetscMatrix( std::declval< const Operator & >().exportMatrix() ) ) >
      {
        using pybind11::operator""_a;

        cls.def_property_readonly( "_backend", [] ( Operator &self ) {
            if (import_petsc4py() != 0)
            {                           \
              std::cout << "ERROR: could not import petsc4py\n";
              throw std::runtime_error("Error during import of petsc4py");
            }
            Mat mat = self.exportMatrix();
            pybind11::handle petsc_mat(PyPetscMat_New(mat));
            return petsc_mat;
          }, pybind11::keep_alive<0,1>());
      }
#endif
#if HAVE_DUNE_ISTL
      template< class B, class A >
      inline static const BCRSMatrix< B, A > &getBCRSMatrix ( const BCRSMatrix< B, A > &matrix ) noexcept
      {
        return matrix;
      }
      template< class Operator, class... options>
      inline static auto addMatrixBackend ( pybind11::class_< Operator, options... > cls, PriorityTag< 3 > )
      -> void_t<decltype( getBCRSMatrix( std::declval< const Operator & >().exportMatrix() ) ) >
      {
        typedef std::decay_t< decltype( getBCRSMatrix( std::declval< const Operator & >().exportMatrix() ) ) > BCRSMatrix;
        if( !pybind11::already_registered< BCRSMatrix >() )
          Python::registerBCRSMatrix< BCRSMatrix >( cls );

        using pybind11::operator""_a;

        cls.def_property_readonly( "_backend", [] ( Operator &self ) {
            return getBCRSMatrix( self.exportMatrix() );
          }, pybind11::keep_alive<0,1>() );
      }
#endif

      template< class AssembledLinearOperator, class... options,
            decltype( std::declval< const AssembledLinearOperator & >().exportMatrix().exportCRS(), 0 ) = 0 >
      inline static void addMatrixBackend ( pybind11::class_< AssembledLinearOperator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def_property_readonly( "_backend", [] ( pybind11::handle self ) {
            auto& mat = self.cast< AssembledLinearOperator &>().exportMatrix();
            auto crs = mat.exportCRS();
            auto &values = std::get<0>(crs);
            auto &inner  = std::get<1>(crs);
            auto &outer  = std::get<2>(crs);
            pybind11::array_t<int> outerIndices(outer.size(),&(outer[0]),self);
            pybind11::array_t<int> innerIndices(inner.size(),&(inner[0]),self);
            pybind11::array_t<double> data(values.size(),&(values[0]),self);
            pybind11::object matrix_type = pybind11::module::import("scipy.sparse").attr("csr_matrix");
            pybind11::object scipy_mat = matrix_type(
                std::tie(data, innerIndices, outerIndices),
                std::make_pair(mat.rows(), mat.cols())
            );
            return scipy_mat;
          }, pybind11::keep_alive<0,1>() );
      }

      // registerLinearOperator
      // ----------------------

      template< class Operator, class... options >
      inline static auto registerLinearOperator ( pybind11::handle scope, pybind11::class_< Operator, options... > cls )
      {
        registerBasicOperator(cls);

#if 0 // linear operators are always assembled
        typedef typename Operator::DomainFunctionType DomainFunction;
        typedef typename Operator::RangeFunctionType RangeFunction;

        // check whether Operator is of type
        // AssembledOperator and thus offers a method matrix.
        static constexpr std::size_t priority =
            std::is_base_of< Dune::Fem::AssembledOperator< DomainFunction, RangeFunction>,
                           Operator > :: value ? 42 : 0;
#endif
        addMatrixBackend(cls, PriorityTag< 42 >());
      }

    } // namespace detail

    // registerOperator
    // ----------------

    template< class Operator, class... options >
    inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
    {
      detail::registerOperator< Operator >( module, cls );
    }

    // registerLinearOperator
    // ----------------------

    template< class Operator, class... options >
    inline static void registerLinearOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
    {
      detail::registerLinearOperator< Operator >( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_OPERATOR_HH
