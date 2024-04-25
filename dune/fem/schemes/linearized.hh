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
        From consturction of the underlying scheme we have
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
    struct LinearScheme : public Scheme::LinearOperatorType
    {
      typedef typename Scheme::LinearOperatorType BaseType;
      typedef Scheme SchemeType;
      typedef typename SchemeType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename SchemeType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename SchemeType::GridPartType GridPartType;
      typedef typename SchemeType::InverseOperatorType LinearInverseOperatorType;
      typedef typename LinearInverseOperatorType::SolverParameterType ParameterType;
      typedef typename SchemeType::ModelType ModelType;

      typedef typename SchemeType::JacobianOperatorType JacobianOperatorType;
      typedef typename SchemeType::DomainFunctionType DomainFunctionType;
      typedef typename SchemeType::RangeFunctionType RangeFunctionType;

      typedef typename SchemeType::DirichletBlockVector DirichletBlockVector;

      struct SolverInfo
      {
        SolverInfo ( bool converged, int linearIterations, int nonlinearIterations )
          : converged( converged ), linearIterations( linearIterations ), nonlinearIterations( nonlinearIterations )
        {}

        bool converged;
        int linearIterations, nonlinearIterations;
        std::vector<double> timing;
      };

      LinearScheme ( SchemeType &scheme,
                     Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : BaseType( "linearized Op", scheme.space(), scheme.space() ),
          scheme_( scheme ),
          inverseOperator_( nullptr ),
          isBound_(false),
          parameter_( std::move( parameter ) ),
          zero_( "zero", scheme.space() )
      {
        zero_.clear();
      }
      LinearScheme ( SchemeType &scheme, const DiscreteFunctionType &ubar,
                     Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : BaseType( "linearized Op", scheme.space(), scheme.space() ),
          scheme_( scheme ),
          inverseOperator_( nullptr ),
          isBound_(false),
          parameter_( std::move( parameter ) ),
          zero_( "zero", scheme.space() )
      {
        zero_.clear();
        scheme.jacobian(ubar,*this);
      }

      /** Note: this sets the error message of the non-existing
       *  non-linear solver and must be here in order to make the
       *  python bindings happy!
       */
      void setErrorMeasure() const {}
      void setConstraints( DomainFunctionType &u ) const
      {
        scheme_.setConstraints(u);
      }
      void setConstraints( const typename DiscreteFunctionType::RangeType &value, DiscreteFunctionType &u ) const
      {
        scheme_.setConstraints(value, u);
      }
      void setConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        scheme_.setConstraints(u, v);
      }
      template < class GridFunctionType,
                 typename = std::enable_if_t< std::is_base_of<Dune::Fem::HasLocalFunction, GridFunctionType>::value > >
      void setConstraints( const GridFunctionType &u, DiscreteFunctionType &v ) const
      {
        scheme_.setConstraints(u, v);
      }
      void subConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        scheme_.subConstraints(u, v);
      }
      void subConstraints( DiscreteFunctionType &v ) const
      {
        scheme_.subConstraints(v);
      }
      void addConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        scheme_.addConstraints(u, v);
      }
      void addConstraints( DiscreteFunctionType &v ) const
      {
        scheme_.addConstraints(v);
      }
      const auto& dirichletBlocks() const
      {
        return scheme_.dirichletBlocks();
      }

      virtual void clear() override
      {
        BaseType::clear();
        if (inverseOperator_)
          inverseOperator_->unbind();
        isBound_ = false;
      }

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const override
      {
        BaseType::operator()(u,w);
      }

      template <class GridFunction>
      void operator() ( const GridFunction &arg, DiscreteFunctionType &dest ) const
      {
        DiscreteFunctionType tmp(dest);
        Fem::interpolate(arg,tmp);
        (*this)(tmp,dest);
      }

      SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution) const
      {
        if (!inverseOperator_)
          inverseOperator_ = std::make_shared<LinearInverseOperatorType>(ParameterType(parameter()));
        if (!isBound_)
        {
          inverseOperator_->bind(*this);
          isBound_ = true;
        }

        DiscreteFunctionType rhs0 = rhs;
        setConstraints(rhs0, solution); // copy the sum to solution for iterative solvers

        (*inverseOperator_)(rhs0, solution );
        return SolverInfo( true, (*inverseOperator_).iterations(), 0);
      }

      /**
       * Solve the system defined by the affine-linear operator
       * without additional rhs, i.e. the rhs is implied by the
       * "affine shift" of the underlying affine linear
       * operator. Dirichlet constraints will be enforced if present
       * in the model.
       */
      SolverInfo solve ( DiscreteFunctionType &solution ) const
      {
        if (!inverseOperator_)
          inverseOperator_ = std::make_shared<LinearInverseOperatorType>(ParameterType(parameter()));
        if (!isBound_)
        {
          inverseOperator_->bind(*this);
          isBound_ = true;
        }
        zero_.clear();
        setConstraints(typename DiscreteFunctionType::RangeType(0), solution);
        (*inverseOperator_)( zero_, solution );
        return SolverInfo( true, (*inverseOperator_).iterations(), 0 );
      }

      bool mark ( double tolerance ) { return scheme_.mark( tolerance ); }
      double estimate ( const DiscreteFunctionType &solution ) { return scheme_.estimate( solution ); }

      const GridPartType &gridPart () const { return scheme_.gridPart(); }
      const DiscreteFunctionSpaceType &space() const { return scheme_.space(); }
      const SchemeType &scheme() const { return scheme_; }

      const ParameterReader& parameter () const { return parameter_; }

    protected:
      SchemeType &scheme_;
      mutable std::shared_ptr<LinearInverseOperatorType> inverseOperator_;
      mutable bool isBound_;
      Dune::Fem::ParameterReader parameter_;
      DiscreteFunctionType zero_;
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
      typedef typename SchemeType::InverseOperatorType LinearInverseOperatorType;
      typedef typename SchemeType::ModelType ModelType;

      typedef typename SchemeType::JacobianOperatorType JacobianOperatorType;
      typedef typename SchemeType::DomainFunctionType DomainFunctionType;
      typedef typename SchemeType::RangeFunctionType RangeFunctionType;

      typedef typename SchemeType::DirichletBlockVector DirichletBlockVector;
      typedef LinearScheme<SchemeType> LinearOperatorType;
      typedef typename LinearOperatorType::SolverInfo SolverInfo;

      LinearizedScheme ( SchemeType &scheme,
                         Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : linOp_( scheme, parameter ),
          rhs_( "affine shift", scheme.space() ),
          ubar_( "ubar", scheme.space() )
      {
        setup();
      }
      LinearizedScheme ( SchemeType &scheme, const DiscreteFunctionType &ubar,
                         Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
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
      void operator() ( const GridFunction &arg, DiscreteFunctionType &dest ) const
      {
        linOp_(arg,dest);
        dest -= rhs_;
      }

      SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution) const
      {
        DiscreteFunctionType sumRhs = rhs;
        sumRhs.axpy(1.0, rhs_);
        // rhs_ = DS[u]u-S[u] and = u-(u-g) = g on boundary so
        // solution = u - DS[u]^{-1}(rhs+S[u]) and = rhs+g on the boundary
        return linOp_.solve( sumRhs, solution);     // don't add constraints again
      }

      /**
       * Solve the system defined by the affine-linear operator
       * without additional rhs, i.e. the rhs is implied by the
       * "affine shift" of the underlying affine linear
       * operator. Dirichlet constraints will be enforced if present
       * in the model.
       */
      SolverInfo solve ( DiscreteFunctionType &solution ) const
      {
        // rhs_ = DS[u]u-S[u] and = u-(u-g) = g on boundary so
        // solution = u - DS[u]^{-1}S[u] and = g on the boundary
        return linOp_.solve( rhs_, solution );
      }

      bool mark ( double tolerance ) { return linOp_.mark( tolerance ); }
      double estimate ( const DiscreteFunctionType &solution ) { return linOp_.estimate( solution ); }

      const GridPartType &gridPart () const { return linOp_.gridPart(); }
      const DiscreteFunctionSpaceType &space() const { return linOp_.space(); }
      const SchemeType &scheme() { return linOp_.scheme(); }
      const ParameterReader& parameter () const { return linOp_.parameter(); }

    protected:
      void setup_(bool isZero=false)
      {
        scheme().jacobian(ubar_, linOp_);

        // compute rhs
        DiscreteFunctionType tmp(ubar_);
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
