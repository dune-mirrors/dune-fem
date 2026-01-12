#ifndef DUNE_FEM_SCHEMES_TUPLEGALERKIN_HH
#define DUNE_FEM_SCHEMES_TUPLEGALERKIN_HH

#include <dune/fem/schemes/galerkin.hh>

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
      : public virtual Operator< typename Op::DomainFunctionType, typename Op::RangeFunctionType >
    {
      typedef typename Op::DomainFunctionType DomainFunctionType;
      typedef typename Op::RangeFunctionType  RangeFunctionType;

      typedef typename RangeFunctionType::GridPartType GridPartType;

      typedef std::tuple< Op,  Args... >  Operators;

      typedef std::make_index_sequence< std::tuple_size< Operators >::value > OperatorIndices;

      static_assert( std::is_same< typename DomainFunctionType::GridPartType, typename RangeFunctionType::GridPartType >::value, "DomainFunction and RangeFunction must be defined on the same grid part." );

      explicit TupleGalerkinOperator ( const Op& op, const Args&... args )
        : ops_( std::make_tuple( op, args... ) )
      {
      }

      void setCommunicate( const bool communicate )
      {
        Hybrid::forEach( OperatorIndices(), [ this, &communicate ] ( auto i ) {
              std::get< i >( this->ops_ ).setCommunicate( communicate );
            });
      }

      void setQuadratureOrders(unsigned int interior, unsigned int surface)
      {
        Hybrid::forEach( OperatorIndices(), [ this, &interior, &surface ] ( auto i ) {
              std::get< i >( this->ops_ ).setQuadratureOrders(interior,surface);
            });
      }

      virtual bool nonlinear() const final override
      {
        bool nonlin = true;
        Hybrid::forEach( OperatorIndices(), [ this, &nonlin ] ( auto i ) {
            nonlin = std::max( nonlin, std::get< i >( this->ops_ ).nonlinear() );
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
      //typedef typename Op::DirichletModelType  DirichletModelType;
      ModelType &model() const { return std::get< 0 >( ops_ ).model(); }

      std::size_t gridSizeInterior () const { return std::get< 0 >( ops_ ).gridSizeInterior(); }

    protected:
      void evaluate( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        w.clear();
        Hybrid::forEach( OperatorIndices(), [ this, &u, &w ] ( auto i ) {
            typedef decltype( std::get< i >( this->ops_ ) ) ThisOp;

            const ThisOp& thisOp = std::get< i >( this->ops_ );
            typedef typename ThisOp::DomainFunction DomFunc;
            typedef typename ThisOp::DomainFunction RanFunc;

            DomFunc d("TupleScheme::eval", thisOp.domainSpace() );
            {
              const auto uit = u.dbegin();
              const auto end = d.dend();
              for( auto it = d.dbegin(); it != end; ++it, ++ uit )
                *it = *uit;
            }
            RanFunc r("TupleScheme::eval", thisOp.rangeSpace() );

            thisOp.evaluate( d, r );

            {
              const auto rit = r.dbegin();
              const auto end = w.dend();
              for( auto it = w.dbegin(); it != end; ++it, ++ rit )
                *it += *rit;
            }
          });
      }

      template < class GridFunction >
      void evaluate( const GridFunction &u, RangeFunctionType &w ) const
      {
        w.clear();
        Hybrid::forEach( OperatorIndices(), [ this, &u, &w ] ( auto i ) {
            typedef decltype( std::get< i >( this->ops_ ) ) ThisOp;

            const ThisOp& thisOp = std::get< i >( this->ops_ );
            typedef typename ThisOp::DomainFunction RanFunc;

            RanFunc r("TupleScheme::eval", thisOp.rangeSpace() );
            thisOp.evaluate( u, r );

            {
              const auto rit = r.dbegin();
              const auto end = w.dend();
              for( auto it = w.dbegin(); it != end; ++it, ++ rit )
                *it += *rit;
            }
          });
      }

      Operators ops_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class Op, class JacobianOperator, class... Args >
    class TupleDifferentiableGalerkinOperator
      : public TupleGalerkinOperator< Op, Args... >,
        public DifferentiableOperator< JacobianOperator >
    {
      typedef TupleGalerkinOperator< Op, Args... >  BaseType;
      typedef typename BaseType::OperatorIndices  OperatorIndices;
      typedef typename BaseType::Operators        Operators;
    public:
      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType    DomainDiscreteFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType     RangeDiscreteFunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;

      explicit TupleDifferentiableGalerkinOperator ( const DomainDiscreteFunctionSpaceType &dSpace,
                                                     const RangeDiscreteFunctionSpaceType &rSpace,
                                                     const Op& op, const Args&... args )
        : BaseType( op, args... )
      {
      }

      void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
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

      void assemble( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
      {
        Hybrid::forEach( OperatorIndices(), [ this, &u, &jOp ] ( auto i ) {
            typedef decltype( std::get< i >( this->ops_ ) ) ThisOp;
            const ThisOp& thisOp = std::get< i >( this->ops_ );

            typedef typename ThisOp::DomainFunction DomFunc;
            typedef typename ThisOp::DomainFunction RanFunc;

            DomFunc d("TupleScheme::assemble", thisOp.domainSpace() );
            {
              const auto uit = u.dbegin();
              const auto end = d.dend();
              for( auto it = d.dbegin(); it != end; ++it, ++ uit )
                *it = *uit;
            }

            thisOp.assemble( d, jOp,
                             /* prepare */  i == 0,
                             /* finalize */ i == std::tuple_size< Operators >::value - 1 );
            });
      }

      template < class GridFunction >
      void assemble( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        Hybrid::forEach( OperatorIndices(), [ this, &u, &jOp ] ( auto i ) {
              std::get< i >( this->ops_ ).assemble( u, jOp,
                                              /* prepare */  i == 0,
                                              /* finalize */ i == std::tuple_size< Operators >::value - 1 );
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

      template< class Op, class LinearInverseOperator, bool addDirichletBC,
                class... OtherOps >
      struct TupleGalerkinSchemeImpl
        : public FemScheme< TupleDifferentiableGalerkinOperator<Op, OtherOps...>, // Operator
                            LinearInverseOperator > // LinearInverseOperator
      {
        typedef FemScheme< TupleDifferentiableGalerkinOperator<Op, OtherOps...>, // Operator
                           LinearInverseOperator >           BaseType;

        typedef typename BaseType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;

        TupleGalerkinSchemeImpl ( const DiscreteFunctionSpaceType &dfSpace,
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
    using VertexCenteredGalerkinScheme = Impl::TupleGalerkinSchemeImpl< DiffOp, InverseOperator, addDirichletBC, AdvOp >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
