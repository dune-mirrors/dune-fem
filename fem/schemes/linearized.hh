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
    // Scheme for
    // L[ubar] + DL[ubar](u-ubar) = 0 in the affine linear case
    // 0 = A ubar - b + Au - A ubar = Au - b = DL[ubar]u - affineShift and
    // so affineShift = b = Aubar - L[ubar]
    template< class Scheme >
    struct LinearizedScheme
    {
      typedef Scheme SchemeType;
      typedef typename SchemeType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename SchemeType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename SchemeType::GridPartType GridPartType;
      typedef typename SchemeType::LinearOperatorType LinearOperatorType;
      typedef typename SchemeType::UsedSolverType::LinearInverseOperatorType LinearInverseOperatorType;
      typedef typename SchemeType::ModelType ModelType;

      LinearizedScheme ( SchemeType &scheme,
                         Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : scheme_( scheme ),
          linearOperator_( "linearized Op", scheme.space(), scheme.space() ),
          linabstol_( 1e-12 ), // parameter.getValue< double >("linabstol", 1e-12 ) ),
          linreduction_( 1e-12 ), // parameter_.getValue< double >("linreduction", 1e-12 ) ),
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
        inverseOperator_ = std::make_shared<LinearInverseOperatorType>(linearOperator_ ,linreduction_, linabstol_ );
        sequence_ = scheme_.space().sequence();
        ubar_.assign(ubar);
      }

      void constraint ( DiscreteFunctionType &u ) const { scheme_.constraint(u); }

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w )
      {
        assert( sequence_ == scheme_.space().sequence() );
        affineShift();
        linearOperator_(u,w);
        w += affineShift_;
      }

      template <class GridFunction>
      void operator() ( const GridFunction &arg, DiscreteFunctionType &dest )
      {
        assert( sequence_ == scheme_.space().sequence() );
        DiscreteFunctionType tmp(dest);
        Fem::interpolate(arg,tmp);
        (*this)(tmp,dest);
      }

      void solve ( DiscreteFunctionType &solution )
      {
        affineShift();
        constraint(solution);
        assert( sequence_ == scheme_.space().sequence() );
        (*inverseOperator_)( affineShift_, solution );
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
      void affineShift()
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
      std::shared_ptr<LinearInverseOperatorType> inverseOperator_;
      DiscreteFunctionType affineShift_;
      Dune::Fem::ParameterReader parameter_;
      DiscreteFunctionType ubar_;
      int sequence_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_LINEARIZED_HH
