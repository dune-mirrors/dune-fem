#ifndef DUNE_FEM_SCHEMES_LINEARIZED_HH
#define DUNE_FEM_SCHEMES_LINEARIZED_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/io/parameter.hh>


namespace Dune
{

  namespace Fem
  {
    // Given non linear operator L[u;(c)] where (c) are some coefficients
    // we define a linear operator
    // A[u;(c)] := L[ubar;(c)] + DL[ubar;(c)](u-ubar)
    //           = DL[ubar;(c)]u + L[ubar;(c)] - DL[ubar;(c)]ubar
    //           = DL[ubar;(c)]u + b
    // Note: DA[ubar,(c)]u = DL[ubar,(c)]u
    // and if L[u];(c)] = Au + c we have
    // A[u;(c)] = Au + Aubar + c - Aubar = Au + c
    //
    // Use in Newton method:
    //   DL[ubar]d + L[ubar] = 0  ; d = -DL^{-1}L[ubar] ; u = u+d
    //   u = u-DL^{-1}L[ubar]
    // A[u] = DL[ubar]u + L[ubar] - DL[ubar]ubar = 0
    //   u = -DL^{-1}(L[ubar]-DL ubar)
    //   u = ubar-DL^{-1}L[ubar]
    template< class Scheme >
    struct LinearizedScheme : public Scheme::LinearOperatorType
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

      LinearizedScheme ( SchemeType &scheme,
                         Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : BaseType( "linearized Op", scheme.space(), scheme.space() ),
          scheme_( scheme ),
          inverseOperator_( nullptr ),
          affineShift_( "affine shift", scheme.space() ),
          isBound_(false),
          parameter_( std::move( parameter ) ),
          ubar_("ubar", scheme.space())
      {
        setup();
      }
      LinearizedScheme ( SchemeType &scheme, const DiscreteFunctionType &ubar,
                         Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : BaseType( "linearized Op", scheme.space(), scheme.space() ),
          scheme_( scheme ),
          inverseOperator_( nullptr ),
          affineShift_( "affine shift", scheme.space() ),
          isBound_(false),
          parameter_( std::move( parameter ) ),
          ubar_("ubar", scheme.space())
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
        setup_();
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
      const auto& dirichletBlocks() const
      {
        return scheme_.dirichletBlocks();
      }

      virtual void clear() override
      {
        BaseType::clear();
        affineShift_.clear();
        inverseOperator_->unbind();
        isBound_ = false;
      }

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
      {
        assert( sequence_ == scheme_.space().sequence() );
        BaseType::operator()(u,w);
        w += affineShift_;
      }

      template <class GridFunction>
      void operator() ( const GridFunction &arg, DiscreteFunctionType &dest ) const
      {
        assert( sequence_ == scheme_.space().sequence() );
        DiscreteFunctionType tmp(dest);
        Fem::interpolate(arg,tmp);
        (*this)(tmp,dest);
      }

      /**
       * additiveConstraints == true:
       * Solve the system with an additional "rhs". That is, if the
       * affine-linear operator is not linear (i.e. already has a
       * "rhs") then the given rhs is added to the already present
       * rhs, including values set for Dirichlet (or other constraint)
       * DOFs. In order to, e.g., solve with Dirichlet zero values you
       * have to install the negative Dirichlet values in the given
       * rhs.
       *
       * additiveConstraints == false
       * Extrawurst in order to be
       * backwards-compatible: Ignore the Dirichlet values contained
       * in rhs. This implies that this variant is not suitable as
       * solver in e.g. Uzawa methods or Newton iterations as it
       * enforces the non-homogeneous Dirichlet values from the model.
       */
      SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution, bool additiveConstraints ) const
      {
        if (!isBound_)
          inverseOperator_->bind(*this);
        DiscreteFunctionType sumRhs = rhs;
        if (!additiveConstraints) {
          setConstraints(typename DiscreteFunctionType::RangeType(0), sumRhs);
        }
        sumRhs.axpy(1.0, affineShift_);
        setConstraints(sumRhs, solution);
        (*inverseOperator_)( sumRhs, solution );
        return SolverInfo( true, (*inverseOperator_).iterations(), 0);
      }

      /**
       * This is the "old" solve with rhs which ignores any values of
       * rhs set for the Constraint DOFs.
       */
      SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution ) const
      {
        return solve(rhs, solution, false);
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
        if (!isBound_)
          inverseOperator_->bind(*this);
        setConstraints(solution);
        (*inverseOperator_)( affineShift_, solution );
        return SolverInfo( true, (*inverseOperator_).iterations(), 0 );
      }

      template< class GridFunction >
      const BaseType &assemble ( const GridFunction &u )
      {
        assert( sequence_ == scheme_.space().sequence() );
        return *this;
      }

      bool mark ( double tolerance ) { return scheme_.mark( tolerance ); }
      double estimate ( const DiscreteFunctionType &solution ) { return scheme_.estimate( solution ); }

      const GridPartType &gridPart () const { return scheme_.gridPart(); }
      const DiscreteFunctionSpaceType &space() const { return scheme_.space(); }

      const ParameterReader& parameter () const
      {
        return parameter_;
      }

    protected:
      void setup_()
      {
        if (!inverseOperator_)
          inverseOperator_ = std::make_shared<LinearInverseOperatorType>(ParameterType(parameter()));
        scheme_.fullOperator().jacobian(ubar_, *this);
        sequence_ = scheme_.space().sequence();
        affineShift();
        inverseOperator_->bind(*this);
        isBound_ = true;
      }

      // this method is called from the 'setup' method and
      // computes 'b' in (with u=ubar):
      // L[w] = DS[u](w-u) + S[u]
      //      = DS[u]w + S[u] - DS[u]u
      //      = DS[u]w + b
      // compute 'b = S[u] - DS[u]u'
      void affineShift() const
      {
        DiscreteFunctionType tmp(ubar_);
        tmp.clear();
        affineShift_.clear();

        // compute DS[u]u
        BaseType::operator()( ubar_, affineShift_ );
        // compute S[u]
        scheme_.fullOperator()( ubar_, tmp );
        affineShift_ -= tmp;
        // ???? affineShift_ *= -1;
        setConstraints(affineShift_);
      }

      SchemeType &scheme_;
      mutable std::shared_ptr<LinearInverseOperatorType> inverseOperator_;
      mutable DiscreteFunctionType affineShift_;
      bool isBound_;
      Dune::Fem::ParameterReader parameter_;
      DiscreteFunctionType ubar_;
      mutable int sequence_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_LINEARIZED_HH
