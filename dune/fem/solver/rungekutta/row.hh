#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_ROW_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_ROW_HH

//- system includes
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-fem includes
#include <dune/fem/solver/rungekutta/basicrow.hh>
#include <dune/fem/solver/rungekutta/butchertable.hh>
#include <dune/fem/solver/rungekutta/timestepcontrol.hh>

namespace DuneODE
{

  // ROWRungeKuttaSolver
  // ------------------------

  /** \brief ROW RungeKutta ODE solver. */
  template< class HelmholtzOperator, class NonlinearSolver, class TimeStepControl = ImplicitRungeKuttaTimeStepControl >
  class ROWRungeKuttaSolver
  : public BasicROWRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl >
  {
    typedef ROWRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl > ThisType;
    typedef BasicROWRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl > BaseType;

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
    ROWRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                          TimeProviderType &timeProvider,
                          int order = 3,
                          const ParametersType &parameters = ParametersType() )
    : BaseType( helmholtzOp, butcherTable( order ), TimeStepControlType( timeProvider, parameters ) )
    {}

  protected:
    static ROWSimpleButcherTable< double > butcherTable ( int order )
    {
      switch( order )
      {
      case 3:
        return row3ButcherTable();
      default:
        DUNE_THROW( NotImplemented, "ROW Runge-Kutta method of order " << order << " not implemented." );
      }
    }
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_ROW_HH
