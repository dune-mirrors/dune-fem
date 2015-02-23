#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_IMPLICIT_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_IMPLICIT_HH

//- system includes
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-fem includes
#include <dune/fem/solver/rungekutta/basicimplicit.hh>
#include <dune/fem/solver/rungekutta/butchertable.hh>
#include <dune/fem/solver/rungekutta/timestepcontrol.hh>

namespace DuneODE
{

  // ImplicitRungeKuttaSolver
  // ------------------------

  /** \brief Implicit RungeKutta ODE solver. */
  template< class HelmholtzOperator, class NonlinearSolver, class TimeStepControl = ImplicitRungeKuttaTimeStepControl >
  class ImplicitRungeKuttaSolver
  : public BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl >
  {
    typedef ImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl > ThisType;
    typedef BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl > BaseType;

  public:
    typedef HelmholtzOperator HelmholtzOperatorType;
    typedef typename BaseType::TimeStepControlType TimeStepControlType;

    typedef typename TimeStepControlType::TimeProviderType TimeProviderType;
    typedef typename TimeStepControlType::ParametersType ParametersType;

    /** \brief constructor
     *
     *  \param[in]  helmholtzOp   Helmholtz operator \f$L\f$
     *  \param[in]  timeProvider  time provider
     *  \param[in]  butcherTable  butcher table to use
     *  \param[in]  traits        additional traits
     */
    ImplicitRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                               TimeProviderType &timeProvider,
                               int order = 1,
                               const ParametersType &parameters = ParametersType() )
    : BaseType( helmholtzOp, butcherTable( order ), TimeStepControlType( timeProvider, parameters ) )
    {}

  protected:
    static SimpleButcherTable< double > butcherTable ( int order )
    {
      switch( order )
      {
      case 1:
        return implicitEulerButcherTable();
      case 2:
        return gauss2ButcherTable();
      case 3:
        return implicit3ButcherTable();
      case 4:
        return implicit34ButcherTable();
      default:
        DUNE_THROW( NotImplemented, "Implicit Runge-Kutta method of order " << order << " not implemented." );
      }
    }
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_IMPLICIT_HH
