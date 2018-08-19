#ifndef DUNE_FEMPY_PY_OPERATOR_HH
#define DUNE_FEMPY_PY_OPERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/pybind11/pybind11.hh>

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

      template< class Operator, class... options, decltype( std::declval< const Operator & >().jacobian( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::JacobianOperatorType >() ), 0 ) = 0 >
      inline static void registerGeneralOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "jacobian", [] ( Operator &self, const GeneralGridFunction< typename Operator::DomainFunctionType > &u, typename Operator::JacobianRangeType &jOp ) { self.jacobian( u, jOp ); }, "u"_a, "jOp"_a );
      }

      template< class Operator, class... options >
      inline static void registerGeneralOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}



      // registerOperatorJacobian
      // ------------------------

      template< class Operator, class... options, decltype( std::declval< const Operator & >().jacobian( std::declval< const typename Operator::DomainFunctionType & >(), std::declval< typename Operator::JacobianOperatorType >() ), 0 ) = 0 >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;

        cls.def( "jacobian", [] ( Operator &self, const typename Operator::DomainFunctionType &u, typename Operator::JacobianRangeType &jOp ) { self.jacobian( u, jOp ); }, "u"_a, "jOp"_a );
      }

      template< class Operator, class... options >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}

      template< class AssembledLinearOperator, class... options>
      inline static void addMatrixBackend ( pybind11::class_< AssembledLinearOperator, options... > cls, PriorityTag< 0 > )
      {
      }

#ifdef PETSC4PY_H // will be set it petsc4py.h was included (so import_petsc4py exists and the python module as well)
      template< class Operator, class... options>
      inline static auto addMatrixBackend ( pybind11::class_< Operator, options... > cls, PriorityTag< 4 > )
      -> void_t<decltype( std::declval<Operator>().petscMatrix() )>
      {
        using pybind11::operator""_a;

        cls.def_property_readonly( "_backend", [] ( Operator &self ) {
            if (import_petsc4py() != 0)
            {                           \
              std::cout << "ERROR: could not import petsc4py\n";
              throw std::runtime_error("Error during import of petsc4py");
            }
            Mat mat = self.petscMatrix();
            pybind11::handle petsc_mat(PyPetscMat_New(mat));
            return petsc_mat;
          });
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
      -> void_t<decltype( getBCRSMatrix( std::declval< const Operator & >().matrix() ) ) >
      {
        typedef std::decay_t< decltype( getBCRSMatrix( std::declval< const Operator & >().matrix() ) ) > BCRSMatrix;
        if( !pybind11::already_registered< BCRSMatrix >() )
          Python::registerBCRSMatrix< BCRSMatrix >( cls );

        using pybind11::operator""_a;

        cls.def_property_readonly( "_backend", [] ( Operator &self ) {
            return getBCRSMatrix( self.matrix() );
          });
      }
#endif

      template< class AssembledLinearOperator, class... options>
      inline static void addMatrixBackend ( pybind11::class_< AssembledLinearOperator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def_property_readonly( "_backend", [] ( AssembledLinearOperator &self ) {
            auto& mat = self.matrix();

            pybind11::array_t<size_t> outerIndices(mat.rows() + 1);

            size_t nnz = 0;
            for(size_t i=0;i<mat.rows();++i)
              nnz += mat.numNonZeros(i);

            pybind11::array_t<size_t> innerIndices(nnz);
            pybind11::array_t<double> data(nnz);

            size_t fill = 0;
            outerIndices.mutable_at(0) = 0;
            for(size_t i=0;i<mat.rows();++i)
            {
              size_t count = i*mat.numNonZeros();
              for(size_t j=0;j<mat.numNonZeros();++j,++count)
              {
                const auto pairIdx = mat.realValue(count);
                if (pairIdx.second < mat.cols())
                {
                  innerIndices.mutable_at(fill) = pairIdx.second;
                  data.mutable_at(fill) = pairIdx.first;
                  ++fill;
                }
                else break;
              }
              outerIndices.mutable_at(i+1) = fill;
            }
            pybind11::object matrix_type = pybind11::module::import("scipy.sparse").attr("csr_matrix");
            pybind11::object scipy_mat = matrix_type(
                std::make_tuple(data, innerIndices, outerIndices),
                std::make_pair(mat.rows(), mat.cols())
            );
            return scipy_mat;
          } );
      }

      template< class GF, class Operator, class... options,
            decltype( std::declval< const Operator & >().jacobian( std::declval< const GF & >(), std::declval< typename Operator::JacobianOperatorType& >() ), 0 ) = 0 >
      inline static void registerOperatorAssemble ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        typedef typename Operator::JacobianOperatorType LinearOperator;

        using pybind11::operator""_a;

        cls.def( "assemble", [] ( const Operator &self, const GF &ubar ) {
            std::unique_ptr<LinearOperator> linOp = std::make_unique<LinearOperator>("tmp", self.domainSpace(), self.rangeSpace());
            self.jacobian( ubar, *linOp );
            return linOp;
          }, "ubar"_a, pybind11::keep_alive<0,1>() );
      }
      template< class GF, class Operator, class... options >
      inline static void registerOperatorAssemble ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {
      }
      template< class Operator, class... options >
      inline static void registerOperatorAssemble ( pybind11::class_< Operator, options... > cls )
      {
        typedef typename Operator::DomainFunctionType DomainFunction;
        typedef typename DomainFunction::RangeType RangeType;
        typedef typename DomainFunction::GridPartType GridPart;

        registerOperatorAssemble< DomainFunction >( cls, PriorityTag< 42 >() );
        registerOperatorAssemble< VirtualizedGridFunction< GridPart, RangeType > >( cls, PriorityTag< 42 >() );
      }

      // registerConstraints
      // -------------------

      template< class Operator, class... options >
      inline static auto registerOperatorConstraints ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval<const Operator>().setConstraints
              ( std::declval<typename Operator::DomainFunctionType&>() ) ) >
      {
        typedef typename Operator::DomainFunctionType DomainFunction;
        typedef typename Operator::RangeFunctionType  RangeFunction;
        cls.def( "setConstraints", [] ( Operator &self, DomainFunction &u) { self.setConstraints( u ); } );
        cls.def( "setConstraints", [] ( Operator &self, const typename DomainFunction::RangeType &value, DomainFunction &u) { self.setConstraints( value, u ); } );
        cls.def( "setConstraints", [] ( Operator &self, const DomainFunction &u, RangeFunction &v) { self.setConstraints( u,v ); } );
        cls.def( "subConstraints", [] ( Operator &self, const DomainFunction &u, RangeFunction &v) { self.subConstraints( u,v ); } );
      }
      template< class Operator, class... options >
      inline static void registerOperatorConstraints ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}
      template< class Operator, class... options >
      inline static void registerOperatorConstraints ( pybind11::class_< Operator, options... > cls )
      {
        registerOperatorConstraints( cls, PriorityTag< 42 >() );
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
      template< class Operator, class... options >
      inline static auto registerJacobianOperator ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      -> void_t< typename Operator::JacobianOperatorType >
      {
        auto clsLinOp = Dune::Python::insertClass< typename Operator::JacobianOperatorType >
              ( cls, "JacobianOperator", Dune::Python::GenerateTypeName(cls,"JacobianOperatorType"));
        if( clsLinOp.second )
        {
          registerBasicOperator(clsLinOp.first);

          typedef typename Operator::DomainFunctionType DomainFunction;
          typedef typename Operator::RangeFunctionType RangeFunction;

          // check whether Operator::JacobianOperatorType is of type
          // AssembledOperator and thus offers a method matrix.
          static constexpr std::size_t priority =
              std::is_base_of< Dune::Fem::AssembledOperator< DomainFunction, RangeFunction>,
                               typename Operator::JacobianOperatorType > :: value ? 42 : 0;
          addMatrixBackend(clsLinOp.first, PriorityTag< priority >());
        }
      }
      template< class Operator, class... options >
      inline static void registerJacobianOperator ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}
      template< class Operator, class... options >
      inline static void registerJacobianOperator ( pybind11::class_< Operator, options... > cls )
      {
        registerJacobianOperator( cls, PriorityTag< 42 >() );
      }


      // registerOperator
      // ----------------

      template< class Operator, class... options >
      inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
      {
        using pybind11::operator""_a;

        registerBasicOperator(cls);

        registerOperatorJacobian( cls, PriorityTag< 42 >() );
        registerGeneralOperatorJacobian( cls, PriorityTag< 42 >() );
        registerOperatorAssemble( cls );
        // registerOperatorModel( cls );
        registerOperatorConstraints( cls );
        registerJacobianOperator( cls );

      }
    } // namespace detail



    // registerOperator
    // ----------------

    template< class Operator, class... options >
    inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
    {
      detail::registerOperator< Operator >( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_OPERATOR_HH
