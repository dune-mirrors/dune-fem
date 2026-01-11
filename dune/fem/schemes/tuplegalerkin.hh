#ifndef DUNE_FEM_SCHEMES_TUPLEGALERKIN_HH
#define DUNE_FEM_SCHEMES_TUPLEGALERKIN_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>
#include <shared_mutex>
#include <vector>
#include <memory>

#include <dune/common/hybridutilities.hh>
#include <dune/common/timer.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/io/parameter/reader.hh>
#include <dune/fem/operator/common/automaticdifferenceoperator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/common/bindguard.hh>

#include <dune/fem/misc/threads/threaditerator.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>
#include <dune/fem/misc/hasboundaryintersection.hh>
#include <dune/fem/misc/griddeclaration.hh>

#include <dune/fem/operator/common/localmatrixcolumn.hh>
#include <dune/fem/operator/common/localcontribution.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/schemes/integrands.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/femscheme.hh>

#include <dune/fem/space/common/capabilities.hh>

// fempy includes
#include <dune/fempy/quadrature/fempyquadratures.hh>

namespace Dune
{

  namespace Fem
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////


    // GalerkinOperator
    // ----------------

    template< class Op, class... Args >
    struct TupleGalerkinOperator
      : public virtual Operator< typename Op::DomainFunction, typename Op::RangeFunction >
    {
      typedef typename Op::DomainFunction DomainFunctionType;
      typedef typename Op::RangeFunction  RangeFunctionType;

      typedef typename RangeFunctionType::GridPartType GridPartType;

      typedef std::tuple< Op,  Args... >  Operators;

      typedef std::make_index_sequence< std::tuple_size< Operators >::value > OperatorIndices;

      static_assert( std::is_same< typename DomainFunctionType::GridPartType, typename RangeFunctionType::GridPartType >::value, "DomainFunction and RangeFunction must be defined on the same grid part." );

      explicit GalerkinOperator ( const Op& op, const Args&... args )
        : ops_( std::make_tuple( op, args ) )
      {
      }

      void setCommunicate( const bool communicate )
      {
        Hybrid::forEach( OperatorIndices(), [ &ops_, &communicate ] ( auto i ) {
              std::get< i >( ops_ ).setCommunicate( communicate );
            });
      }

      void setQuadratureOrders(unsigned int interior, unsigned int surface)
      {
        Hybrid::forEach( OperatorIndices(), [ &ops_, &interior, &surface ] ( auto i ) {
              std::get< i >( ops_ ).setQuadratureOrders(interior,surface);
            });
      }

      virtual bool nonlinear() const final override
      {
        bool nonlin = true;
        Hybrid::forEach( OperatorIndices(), [ &ops_, &nonlin ] ( auto i ) {
            nonlin = std::max( nonlin, std::get< i >( ops_ ).nonlinear() );
            });
        return nonlin;
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const final override
      {
        evaluate( u, w );
      }

      template< class GridFunction >
      void operator() ( const GridFunction &u, RangeFunctionType &w ) const
      {
        evaluate( u, w );
      }

      const GridPartType &gridPart () const { return std::get< 0 >( ops_ ).gridPart(); }

      typedef typename Op::ModelType           ModelType;
      typedef typename Op::DirichletModelType  DirichletModelType;
      ModelType &model() const { return std::get< 0 >( ops_ ).model(); }

      std::size_t gridSizeInterior () const { return std::get< 0 >( ops_ ).gridSizeInterior(); }

    protected:
      template < class GridFunction >
      void evaluate( const GridFunction &u, RangeFunctionType &w ) const
      {
        Hybrid::forEach( OperatorIndices(), [ &ops_, &nonlin ] ( auto i ) {
            nonlin = std::max( nonlin, std::get< i >( ops_ ).nonlinear() );
            });
        iterators_.update();
      }

      Operators ops_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class Op, class, JacobianOperator, class... Args, >
    class TupleDifferentiableGalerkinOperator
      : public TupleGalerkinOperator< Op, Args... >,
        public DifferentiableOperator< JacobianOperator >
    {
      typedef TupleGalerkinOperator< Op, Args... >  BaseType;
      typedef typename BaseType::OperatorIndices  OperatorIndices;
    public:
      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType    DomainDiscreteFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType     RangeDiscreteFunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;

      explicit DifferentiableGalerkinOperator ( const DomainDiscreteFunctionSpaceType &dSpace,
                                                const RangeDiscreteFunctionSpaceType &rSpace,
                                                const Op& op, const Args&... args )
        : BaseType( op, args... )
      {
      }

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const final override
      {
        assemble( u, jOp );
      }

      template< class GridFunction >
      void jacobian ( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        assemble( u, jOp );
      }

      const DomainDiscreteFunctionSpaceType& domainSpace() const
      {
        return std::get< 0 >( ops_ ).domainSpace();
      }
      const RangeDiscreteFunctionSpaceType& rangeSpace() const
      {
        return std::get< 0 >( ops_ ).rangeSpace();
      }

      using BaseType::nonlinear;

    protected:
      using BaseType::ops_;

      template < class GridFunction >
      void assemble( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        Hybrid::forEach( OperatorIndices(), [ &ops_, &u, &jOp ] ( auto i ) {
              std::get< i >( ops_ ).assemble( u, jOp,
                                              /* prepare */  i == 0,
                                              /* finalize */ i == std::tuple_size< Operators > - 1 );
            });
      }
    };


    namespace Impl
    {

      /*
      // VertexCenteredGalerkinSchemeImpl
      // --------------------------------
      template< class DiffOp, class AdvOp, class LinearOperator, bool addDirichletBC >
      struct GalerkinSchemeTraits
      {
        template <class O, bool addDBC>
        struct DirichletBlockSelector { using type = void; };
        template <class O>
        struct DirichletBlockSelector<O,true> { using type = typename O::DirichletBlockVector; };

        using DifferentiableOperatorType = std::conditional_t< addDirichletBC,
           DirichletWrapperOperator< DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator >>,
           DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator > >;
        using DirichletBlockVector = typename DirichletBlockSelector<
                 DirichletWrapperOperator<
                    DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator >>,
                 addDirichletBC>::type;

        typedef DifferentiableOperatorType type;
      };
      */

      template< class Op, class LinearOperator, class LinearInverseOperator, bool addDirichletBC,
                class... OtherOps,
                template <class,class...> class DifferentiableGalerkinOperatorImpl = TupleDifferentiableGalerkinOperator >
      struct TupleGalerkinSchemeImpl
        : public FemScheme< TupleGalerkinOperator<Op, OtherOps...>, // Operator
                            LinearInverseOperator > // LinearInverseOperator
      {
        typedef FemScheme< TupleGalerkinOperator<Op, OtherOps...>, // Operator
                           LinearInverseOperator >           BaseType;

        typedef typename BaseType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;

        GalerkinSchemeImpl ( const DiscreteFunctionSpaceType &dfSpace,
                             const Op& op, const OtherOps&... ops,
                             const ParameterReader& parameter = Parameter::container() )
          : BaseType(dfSpace,
                     parameter,
                     op, ops...)
        {}
      };

    } // end namespace Impl

    // VertexCenteredGalerkinScheme
    // ----------------------------

    template<class DiffOp, class AdvOp, class InverseOperator, bool addDirichletBC >
    using VertexCenteredGalerkinScheme = Impl::TupleGalerkinSchemeImpl< DiffOp, InverseOperator, addDirichletBC, AdvOp, DifferentiableGalerkinOperator >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
