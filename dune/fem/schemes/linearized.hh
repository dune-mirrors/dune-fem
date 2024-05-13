#ifndef DUNE_FEM_SCHEMES_LINEARIZED_HH
#define DUNE_FEM_SCHEMES_LINEARIZED_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/schemes/femscheme.hh>


namespace Dune
{

  namespace Fem
  {
    /*
      Idea is to setup the first two terms in the Taylor expansion:
        L[w] = S[u] + DS[u](w-u)
             = S[u] - DS[u]u + DS[u]w
             = DS[u]w + b

        So b = S[u] - DS[u]u
        and L[w] = DS[u]w + b

        The solve method with rhs=f computs the solution w to
        L[w] = f <=> DS[u]w = f-b

        With Dirichlet BCs:
        From construction of the underlying scheme we have
        S[u] = u - g and DS[u]w = w on the boundary.

        Therefore on the boundary we have
        b = S[u] - DS[u]u = u - g - u = -g
        and L[w] = DS[u]w + b = w - g
        and solving L[w] = f is the same as DS[u]w = f-b or w = f-b = f+g on the boundary

        Note: we store -b in the LinearizedScheme, the LinearScheme
              represents the above but with b=0
    */

    //////////////////////////////////////////////////////////////////////////////////

    // This class implements the Jacobian part of a scheme, i.e., L[w] = DS[u]w
    // without the rhs 'b'.
    template< class Scheme >
    struct LinearScheme : public Scheme::LinearOperatorType,
                          public FemScheme< Scheme, typename Scheme::LinearInverseOperatorType, typename Scheme::LinearInverseOperatorType >
    {
      typedef Scheme SchemeType;
      typedef typename SchemeType::LinearInverseOperatorType LinearInverseOperatorType;

      // base types
      typedef typename SchemeType::LinearOperatorType BaseType;
      typedef FemScheme< SchemeType, LinearInverseOperatorType, LinearInverseOperatorType > FSBaseType;

      typedef typename SchemeType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename SchemeType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename SchemeType::GridPartType GridPartType;
      typedef typename LinearInverseOperatorType::SolverParameterType ParameterType;
      typedef typename SchemeType::ModelType             ModelType;

      typedef typename SchemeType::JacobianOperatorType  JacobianOperatorType;
      typedef typename SchemeType::DomainFunctionType    DomainFunctionType;
      typedef typename SchemeType::RangeFunctionType     RangeFunctionType;

      typedef typename SchemeType::DirichletBlockVector  DirichletBlockVector;

      typedef Dune::Fem::PreconditionerFunctionWrapper<RangeFunctionType,DomainFunctionType >  PreconditionerFunctionWrapperType;
      // std::function to represents the Python function passed as potential preconditioner
      typedef typename PreconditionerFunctionWrapperType::PreconditionerFunctionType  PreconditionerFunctionType;

      using FSBaseType :: fullOperator;
      using FSBaseType :: setConstraints;
      using FSBaseType :: space;

      // from femscheme.hh
      typedef typename LinearInverseOperatorType::SolverInfoType  SolverInfoType ;

      //! \todo @docme
      LinearScheme ( SchemeType &scheme,
                     Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : BaseType( "linearized Op", scheme.space(), scheme.space() ),
          FSBaseType( scheme, parameter ),
          isBound_(false),
          parameter_( std::move( parameter ) ),
          tmp_( "LS::tmp", scheme.space() )
      {
      }

      //! \todo @docme
      LinearScheme ( SchemeType &scheme, const DiscreteFunctionType &ubar,
                     Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : BaseType( "linearized Op", scheme.space(), scheme.space() ),
          FSBaseType( scheme, parameter ),
          isBound_(false),
          parameter_( std::move( parameter ) ),
          tmp_( "LS::tmp", scheme.space() )
      {
        fullOperator().jacobian(ubar,*this);
      }

      /** Note: this sets the error message of the non-existing
       *  non-linear solver and must be here in order to make the
       *  python bindings happy!
       */
      void setErrorMeasure() const {}

      using BaseType::clear;
      virtual void clear() override
      {
        BaseType::clear();
        invOp_.unbind();
        isBound_ = false;
      }

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const override
      {
        // use linear op (instead of fullOperator)
        BaseType::operator()(u,w);
      }

      template <class GridFunction>
      void operator() ( const GridFunction &arg, DiscreteFunctionType &dest ) const
      {
        Fem::interpolate(arg, tmp_);
        // apply operator (i.e. matvec)
        (*this)(tmp_, dest);
      }

      SolverInfoType solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution) const
      {
        if (!isBound_)
        {
          invOp_.bind(*this);
          isBound_ = true;
        }

        // copy Dirichlet constraints from rhs to solution
        setConstraints(rhs, solution);
        invOp_(rhs, solution );
        return invOp_.info();
      }

      SolverInfoType solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution, const PreconditionerFunctionType& p) const
      {
        if (isBound_)
          invOp_.unbind();

        PreconditionerFunctionWrapperType pre( p );
        invOp_.bind(*this, pre);
        isBound_ = true;
        auto info = solve( rhs, solution );
        invOp_.unbind();
        isBound_ = false;
        return info;
      }

      /**
       * Solve the system defined by the affine-linear operator
       * without additional rhs, i.e. the rhs is implied by the
       * "affine shift" of the underlying affine linear
       * operator. Dirichlet constraints will be enforced if present
       * in the model.
       */
      SolverInfoType solve ( DiscreteFunctionType &solution ) const
      {
        if (!isBound_)
        {
          invOp_.bind(*this);
          isBound_ = true;
        }
        DiscreteFunctionType& zero = tmp_;
        zero.clear();
        setConstraints(typename DiscreteFunctionType::RangeType(0), solution);
        invOp_( zero, solution );
        return invOp_.info();
      }

      const SchemeType &scheme() const { return fullOperator(); }
      const ParameterReader& parameter () const { return parameter_; }

      DiscreteFunctionType& temporaryData() const { return tmp_; }
    protected:
      using FSBaseType :: invOp_;
      mutable bool isBound_;
      Dune::Fem::ParameterReader parameter_;
      mutable DiscreteFunctionType tmp_;
    };

    // This class stores a 'LinearScheme' and defines 'setup' methods that
    // compute the rhs 'b'. All methods are then
    template< class Scheme >
    struct LinearizedScheme
    : public Dune::Fem::Operator<typename Scheme::DomainFunctionType, typename Scheme::RangeFunctionType>
    {
      typedef Scheme SchemeType;
      typedef typename SchemeType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename SchemeType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename SchemeType::GridPartType GridPartType;
      typedef typename SchemeType::LinearInverseOperatorType LinearInverseOperatorType;
      typedef typename SchemeType::ModelType ModelType;

      typedef typename SchemeType::JacobianOperatorType JacobianOperatorType;
      typedef typename SchemeType::DomainFunctionType DomainFunctionType;
      typedef typename SchemeType::RangeFunctionType RangeFunctionType;

      typedef typename SchemeType::DirichletBlockVector DirichletBlockVector;
      typedef LinearScheme<SchemeType> LinearOperatorType;
      typedef typename LinearOperatorType::SolverInfoType SolverInfoType;

      typedef Dune::Fem::PreconditionerFunctionWrapper<
          typename LinearOperatorType::RangeFunctionType,
          typename LinearOperatorType::DomainFunctionType >  PreconditionerFunctionWrapperType;
      // std::function to represents the Python function passed as potential preconditioner
      typedef typename PreconditionerFunctionWrapperType::PreconditionerFunctionType  PreconditionerFunctionType;

      LinearizedScheme ( SchemeType &scheme,
                         const Dune::Fem::ParameterReader& parameter = Dune::Fem::Parameter::container() )
        : linOp_( scheme, parameter ),
          rhs_( "affine shift", scheme.space() ),
          ubar_( "ubar", scheme.space() )
      {
        setup();
      }
      LinearizedScheme ( SchemeType &scheme, const DiscreteFunctionType &ubar,
                         const Dune::Fem::ParameterReader& parameter = Dune::Fem::Parameter::container() )
        : linOp_( scheme, ubar, parameter ),
          rhs_( "affine shift", scheme.space() ),
          ubar_( "ubar", scheme.space() )
      {
        setup(ubar);
      }

      void setup(const DiscreteFunctionType &ubar)
      {
        ubar_.assign(ubar);
        setup_();
      }

      template <class GridFunction>
      void setup(const GridFunction &ubar)
      {
        Fem::interpolate(ubar, ubar_);
        setup_();
      }
      void setup()
      {
        ubar_.clear();
        setup_(true); // provide info that `ubar=0`
      }

      /** Note: this sets the error message of the non-existing
       *  non-linear solver and must be here in order to make the
       *  python bindings happy!
       */
      void setErrorMeasure() const {}
      void setConstraints( DomainFunctionType &u ) const
      {
        linOp_.setConstraints(u);
      }
      void setConstraints( const typename DiscreteFunctionType::RangeType &value, DiscreteFunctionType &u ) const
      {
        linOp_.setConstraints(value, u);
      }
      void setConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        linOp_.setConstraints(u, v);
      }
      template < class GridFunctionType,
                 typename = std::enable_if_t< std::is_base_of<Dune::Fem::HasLocalFunction, GridFunctionType>::value > >
      void setConstraints( const GridFunctionType &u, DiscreteFunctionType &v ) const
      {
        linOp_.setConstraints(u, v);
      }
      void subConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        linOp_.subConstraints(u, v);
      }
      void subConstraints( DiscreteFunctionType &v ) const
      {
        linOp_.subConstraints(v);
      }
      void addConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        linOp_.addConstraints(u, v);
      }
      void addConstraints( DiscreteFunctionType &v ) const
      {
        linOp_.addConstraints(v);
      }
      const auto& dirichletBlocks() const
      {
        return linOp_.dirichletBlocks();
      }

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const override
      {
        linOp_(u,w);
        w -= rhs_;
      }
      template <class GridFunction>
      void operator() ( const GridFunction &arg, DiscreteFunctionType &w ) const
      {
        linOp_(arg, w);
        w -= rhs_;
      }

      SolverInfoType solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution) const
      {
        return _solve(rhs, solution, nullptr );
      }

      SolverInfoType solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution, const PreconditionerFunctionType& p) const
      {
        // TODO: Pass preconditioning to solver !
        return _solve( rhs, solution, &p );
      }

      /**
       * Solve the system defined by the affine-linear operator
       * without additional rhs, i.e. the rhs is implied by the
       * "affine shift" of the underlying affine linear
       * operator. Dirichlet constraints will be enforced if present
       * in the model.
       */
      SolverInfoType solve ( DiscreteFunctionType &solution ) const
      {
        // rhs_ = DS[u]u-S[u] and = u-(u-g) = g on boundary so
        // solution = u - DS[u]^{-1}S[u] and = g on the boundary
        return linOp_.solve( rhs_, solution );
      }

      const GridPartType &gridPart () const { return linOp_.gridPart(); }
      const DiscreteFunctionSpaceType &space() const { return linOp_.space(); }
      const SchemeType &scheme() { return linOp_.scheme(); }
      const ParameterReader& parameter () const { return linOp_.parameter(); }

    protected:
      SolverInfoType _solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution, const PreconditionerFunctionType* p) const
      {
        DiscreteFunctionType& sumRhs = linOp_.temporaryData();
        sumRhs.assign(rhs);
        sumRhs += rhs_;

        // rhs_ = DS[u]u-S[u] and = u-(u-g) = g on boundary so
        // solution = u - DS[u]^{-1}(rhs+S[u]) and = rhs+g on the boundary
        if( p )
          return linOp_.solve( sumRhs, solution, *p); // don't add constraints again
        else
          return linOp_.solve( sumRhs, solution);     // don't add constraints again
      }

      void setup_(bool isZero=false)
      {
        scheme().jacobian(ubar_, linOp_);

        // compute rhs
        DiscreteFunctionType& tmp = linOp_.temporaryData();

        // compute tmp = S[u] (tmp = u-g on boundary)
        tmp.clear();
        scheme()( ubar_, tmp );

        // compute rhs_ = DS[u]u (rhs_ = u on boundary)
        rhs_.clear();
        if (!isZero)
          linOp_( ubar_, rhs_ );

        // compute rhs_ - tmp = DS[u]u-S[u] and u-(u-g)=g on boundary
        rhs_ -= tmp;
      }

      LinearOperatorType linOp_;
      DiscreteFunctionType rhs_;
      DiscreteFunctionType ubar_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_LINEARIZED_HH
