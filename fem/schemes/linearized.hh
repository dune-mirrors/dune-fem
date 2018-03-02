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
    struct LinearizedScheme
    {
      typedef Scheme SchemeType;
      typedef typename SchemeType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename SchemeType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename SchemeType::GridPartType GridPartType;
      typedef typename SchemeType::LinearOperatorType LinearOperatorType;
      typedef typename SchemeType::InverseOperatorType LinearInverseOperatorType;
      typedef typename SchemeType::ModelType ModelType;
      struct SolverInfo
      {
        SolverInfo ( bool converged, int linearIterations, int nonlinearIterations )
          : converged( converged ), linearIterations( linearIterations ), nonlinearIterations( nonlinearIterations )
        {}

        bool converged;
        int linearIterations, nonlinearIterations;
      };

      LinearizedScheme ( SchemeType &scheme,
                         Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : scheme_( scheme ),
          linearOperator_( "linearized Op", scheme.space(), scheme.space() ),
          linabstol_( 1e-12 ), // parameter.getValue< double >("linabstol", 1e-12 ) ),
          linreduction_( 1e-12 ), // parameter_.getValue< double >("linreduction", 1e-12 ) ),
          maxIter_( 1000 ),
          inverseOperator_( nullptr ),
          affineShift_( "affine shift", scheme.space() ),
          parameter_( std::move( parameter ) ),
          ubar_("ubar", scheme.space())
      {
        ubar_.clear();
        setup(ubar_);
      }
      LinearizedScheme ( SchemeType &scheme, const DiscreteFunctionType &ubar,
                         Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : scheme_( scheme ),
          linearOperator_( "linearized Op", scheme.space(), scheme.space() ),
          linabstol_( 1e-12 ), // parameter.getValue< double >("linabstol", 1e-12 ) ),
          linreduction_( 1e-12 ), // parameter_.getValue< double >("linreduction", 1e-12 ) ),
          maxIter_( 1000 ),
          inverseOperator_( nullptr ),
          affineShift_( "affine shift", scheme.space() ),
          parameter_( std::move( parameter ) ),
          ubar_("ubar", scheme.space())
      {
        setup(ubar);
      }

      void setup(const DiscreteFunctionType &ubar)
      {
        scheme_.fullOperator().jacobian(ubar, linearOperator_);
        inverseOperator_ = std::make_shared<LinearInverseOperatorType>(linreduction_, linabstol_, maxIter_ );
        inverseOperator_->bind(linearOperator_);
        sequence_ = scheme_.space().sequence();
        ubar_.assign(ubar);
      }

      void constraint ( DiscreteFunctionType &u ) const { scheme_.constraint(u); }

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
      {
        assert( sequence_ == scheme_.space().sequence() );
        affineShift();
        linearOperator_(u,w);
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

      SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution ) const
      {
        int oldCount = (*inverseOperator_).iterations();
        constraint(solution);
        assert( sequence_ == scheme_.space().sequence() );
        (*inverseOperator_)( rhs, solution );
        return SolverInfo( true, (*inverseOperator_).iterations()-oldCount, 1 );
      }
      SolverInfo solve ( DiscreteFunctionType &solution ) const
      {
        affineShift();
        return solve(affineShift_,solution);
      }

      template< class GridFunction >
      const LinearOperatorType &assemble ( const GridFunction &u )
      {
        assert( sequence_ == scheme_.space().sequence() );
        return linearOperator_;
      }

      bool mark ( double tolerance ) { return scheme_.mark( tolerance ); }
      double estimate ( const DiscreteFunctionType &solution ) { return scheme_.estimate( solution ); }

      const GridPartType &gridPart () const { return scheme_.gridPart(); }
      const DiscreteFunctionSpaceType &space() const { return scheme_.space(); }

    protected:
      void affineShift() const
      {
        DiscreteFunctionType tmp(ubar_);
        tmp.clear();

        affineShift_.clear();
        linearOperator_( ubar_, affineShift_ );
        scheme_.fullOperator()( ubar_, tmp );
        affineShift_ -= tmp;
        constraint(affineShift_);
      }

      SchemeType &scheme_;
      LinearOperatorType linearOperator_;
      double linabstol_;
      double linreduction_;
      int maxIter_;
      std::shared_ptr<LinearInverseOperatorType> inverseOperator_;
      mutable DiscreteFunctionType affineShift_;
      Dune::Fem::ParameterReader parameter_;
      DiscreteFunctionType ubar_;
      int sequence_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_LINEARIZED_HH
