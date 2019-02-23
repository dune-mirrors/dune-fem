#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_ROW_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_ROW_HH

//- system includes
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-fem includes
#include <dune/fem/io/parameter.hh>
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

    typedef typename TimeStepControlType::TimeProviderType   TimeProviderType;
    typedef typename TimeStepControlType::ParameterType      TimeStepControlParameterType;
    typedef typename BaseType::NonlinearSolverParameterType  NonlinearSolverParameterType;
    typedef NonlinearSolverParameterType  ParameterType;

    /** \brief constructor
     *
     *  \param[in]  helmholtzOp   Helmholtz operator \f$L\f$
     *  \param[in]  timeProvider  time provider
     *  \param[in]  order         order of butcher table to use
     *  \param[in]  tscParam      parameters for implicit time step control
     *  \param[in]  nlsParam      parameters for non linear solver control
     */
    ROWRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                          TimeProviderType &timeProvider, int order,
                          const TimeStepControlParameterType& tscParam,
                          const NonlinearSolverParameterType& parameter )
    : BaseType( helmholtzOp, timeProvider, butcherTable( order ), TimeStepControlType( timeProvider, tscParam ), parameter )
    {}

    ROWRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                          TimeProviderType &timeProvider, int order,
                          const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType( helmholtzOp, timeProvider, butcherTable( order ), TimeStepControlType( timeProvider, parameter ),
                NonlinearSolverParameterType( parameter ) )
    {}

    /** \brief constructor
     *
     *  \param[in]  helmholtzOp   Helmholtz operator \f$L\f$
     *  \param[in]  timeProvider  time provider
     *  \param[in]  tscParam      parameters for implicit time step control
     *  \param[in]  nlsParam      parameters for non linear solver control
     */
    ROWRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                          TimeProviderType &timeProvider,
                          const TimeStepControlParameterType& tscParam,
                          const NonlinearSolverParameterType& parameter )
    : BaseType( helmholtzOp, timeProvider, butcherTable( 3 ), TimeStepControlType( timeProvider, tscParam ), parameter )
    {}

    ROWRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                          TimeProviderType &timeProvider,
                          const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType( helmholtzOp, timeProvider, butcherTable( 3 ), TimeStepControlType( timeProvider, parameter ),
        ParameterType( parameter ), NonlinearSolverParameterType( parameter ) )
    {}

  protected:
    static ROWSimpleButcherTable< double > butcherTable ( int order )
    {
      switch( order )
      {
      case 2:
        return row2ButcherTable();
      case 3:
        return row3ButcherTable();
      default:
        DUNE_THROW( NotImplemented, "ROW Runge-Kutta method of order " << order << " not implemented." );
      }
    }
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_ROW_HH
