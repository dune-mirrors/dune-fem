#ifndef DUNE_FEM_TIMEDEPENDENTFUNCTION_HH
#define DUNE_FEM_TIMEDEPENDENTFUNCTION_HH

#include <dune/common/deprecated.hh>

#include <dune/fem/function/common/instationary.hh>

namespace Dune
{

  namespace Fem
  {

    // TimeDependentFunction
    // ---------------------

    template< class Function >
    class TimeDependentFunction
      : public InstationaryFunction< Function, __InstationaryFunction::HoldReference >
    {
      typedef InstationaryFunction< Function, __InstationaryFunction::HoldReference > BaseType;

    public:
      TimeDependentFunction ( const Function &function, double time )
        : BaseType( function, time )
      {}

      /** \copydoc Dune::Fem::Function::evaluate */
      void evaluate ( const typename BaseType::DomainType &x,
                      typename BaseType::RangeType &value ) const
      DUNE_DEPRECATED_MSG( "Dune::Fem::TimeDependentFunction has been deprecated, use Dune::Fem::InstationaryFunction instead" )
      {
        BaseType::evaluate( x, value );
      }

      /** \copydoc Dune::Fem::Function::jacobian */
      void jacobian ( const typename BaseType::DomainType &x,
                      typename BaseType::JacobianRangeType &jacobian ) const
      DUNE_DEPRECATED_MSG( "Dune::Fem::TimeDependentFunction has been deprecated, use Dune::Fem::InstationaryFunction instead" )
      {
        BaseType::jacobian( x, jacobian );
      }

      /** \copydoc Dune::Fem::Function::hessian */
      void hessian ( const typename BaseType::DomainType &x,
                     typename BaseType::HessianRangeType &hessian ) const
      DUNE_DEPRECATED_MSG( "Dune::Fem::TimeDependentFunction has been deprecated, use Dune::Fem::InstationaryFunction instead" )
      {
        BaseType::hessian( x, hessian );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_HH
