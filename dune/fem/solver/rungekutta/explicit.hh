#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_EXPLICIT_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_EXPLICIT_HH

//- system includes
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

//- Dune includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/odesolverinterface.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/solver/rungekutta/butchertable.hh>

namespace DuneODE
{

  using namespace Dune;
  using namespace Fem;
  using namespace std;

  /** @addtogroup ODESolver

      A variety of runge kutta have so far been implemented based
      on the \c pardg package.
      These include explicit ssp methods upto order 4
      through Dune::ExplicitRungeKuttaSolver and Dune::ExplicitTimeStepper
      (where the second version uses \c pardg).
      Furthermore implicit methods and IMEX schemes up to order
      4 can be used through Dune::ImplicitOdeSolver
      and Dune::SemiImplicitOdeSolver.
      Each of these classes require an operator defining the right hand
      side for its construction; the template argument is a
      Dune::DiscreteFunctionSpace and the operator has to be derived
      from the Dune::SpaceOperatorInterface<DiscreteFunctionSpace>.
      Further arguments in the constructor are an instance of the
      Dune::Timeprovider and the requested order of the scheme which can
      be chosen during runtime.

      The management of the simulation time is performed by an instance of the
      Dune::TimeProvider class. The time at which the space operator is
      to be evaluated is passed to the operator in the \c setTime method on
      the Dune::SpaceOperatorInterface and it is assumed that
      through the method \c timeStepEstimate on the
      Dune::SpaceOperatorInterface an upper estimate for the time
      step is provided satifying at the least the stability
      restriction of the forward euler method. This estimate is then
      multiplied by an additional factor depending on the ode solver used
      and then passed to the instance of the Dune::TimeProvider.

      \remarks
      The interface for ODE solvers is defined by the class
      OdeSolverInterface. The interface for discretization operators
      working with the OdeSolvers is described by the class SpaceOperatorInterface.
   @{
   **/


  /** \brief Exlicit RungeKutta ODE solver. */
  template<class DestinationImp>
  class ExplicitRungeKuttaSolver
  : public OdeSolverInterface<DestinationImp>
  {
  public:
    typedef DestinationImp DestinationType;
    typedef SpaceOperatorInterface<DestinationImp> OperatorType;
    typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;

    typedef typename OdeSolverInterface<DestinationImp> :: MonitorType MonitorType ;

    using OdeSolverInterface<DestinationImp> :: solve ;
  protected:
    SimpleButcherTable< double > defaultButcherTables( const int order ) const
    {
      switch( order )
      {
        case 1: return explicitEulerButcherTable();
        case 2: return tvd2ButcherTable();
        case 3: return tvd3ButcherTable();
        case 4: return rk4ButcherTable();
        case 5:
        case 6: return expl6ButcherTable();

        default:
               std::cerr<< "Warning: ExplicitRungeKutta of order "<< order << " not implemented, using order 6!" << std::endl;
               return expl6ButcherTable();
      }
    }

  public:
    /** \brief constructor
      \param[in] op Operator \f$L\f$
      \param[in] tp TimeProvider
      \param[in] bt Butcher table defining the Runge-Kutta scheme
      \param[in] verbose verbosity
    */
    ExplicitRungeKuttaSolver(OperatorType& op,
                             TimeProviderBase& tp,
                             const SimpleButcherTable< double >& butcherTable,
                             bool verbose )
      : A_( butcherTable.A() ),
        b_( butcherTable.b() ),
        c_( butcherTable.c() ),
        Upd(),
        op_(op),
        tp_(tp),
        ord_( butcherTable.order() ),
        stages_( butcherTable.stages() ),
        initialized_(false)
    {
      assert(ord_>0);

      // create update memory
      for (int i=0; i<stages_; ++i)
      {
        Upd.emplace_back( new DestinationType("URK",op_.space()) );
      }
      Upd.emplace_back(new DestinationType("Ustep",op_.space()) );

    }

    /** \brief constructor
      \param[in] op Operator \f$L\f$
      \param[in] tp TimeProvider
      \param[in] pord polynomial order
      \param[in] verbose verbosity
    */
    ExplicitRungeKuttaSolver(OperatorType& op,
                             TimeProviderBase& tp,
                             const int pord,
                             bool verbose )
      : ExplicitRungeKuttaSolver( op, tp, defaultButcherTables( pord ), verbose )
    {
    }

    ExplicitRungeKuttaSolver(OperatorType& op,
                             TimeProviderBase& tp,
                             const int pord,
                             const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : ExplicitRungeKuttaSolver( op, tp, pord, true )
    {}

    //! apply operator once to get dt estimate
    void initialize(const DestinationType& U0)
    {
      if( ! initialized_ )
      {
        // Compute Steps
        op_(U0, *(Upd[0]));
        initialized_ = true;

        // provide operators time step estimate
        tp_.provideTimeStepEstimate( op_.timeStepEstimate() );
      }
    }

    //! solve the system
    void solve(DestinationType& U0, MonitorType& monitor )
    {
      // no information here
      monitor.reset();

      // initialize
      if( ! initialized_ )
      {
        DUNE_THROW(InvalidStateException,"ExplicitRungeKuttaSolver wasn't initialized before first call!");
      }

      // get cfl * timeStepEstimate
      const double dt = tp_.deltaT();
      // get time
      const double t = tp_.time();

      // set new time
      op_.setTime( t );

      // Compute Steps
      op_(U0, *(Upd[0]));

      // provide operators time step estimate
      tp_.provideTimeStepEstimate( op_.timeStepEstimate() );

      for (int i=1; i<stages_; ++i)
      {
        (Upd[ord_])->assign(U0);
        for (int j=0; j<i ; ++j)
        {
          (Upd[ord_])->axpy((A_[i][j]*dt), *(Upd[j]));
        }

        // set new time
        op_.setTime( t + c_[i]*dt );

        // apply operator
        op_( *(Upd[ord_]), *(Upd[i]) );

        // provide operators time step estimate
        tp_.provideTimeStepEstimate( op_.timeStepEstimate() );
      }

      // Perform Update
      for (int j=0; j<stages_; ++j)
      {
        U0.axpy((b_[j]*dt), *(Upd[j]));
      }
    }

    void description(std::ostream& out) const
    {
      out << "ExplRungeKutta, steps: " << ord_
          //<< ", cfl: " << this->tp_.factor()
          << "\\\\" <<std::endl;
    }

  protected:
    // Butcher table A,b,c
    Dune::DynamicMatrix< double > A_;
    Dune::DynamicVector< double > b_;
    Dune::DynamicVector< double > c_;

    // stages of Runge-Kutta solver
    std::vector< std::unique_ptr< DestinationType > > Upd;

    // operator to solve for
    OperatorType& op_;
    // time provider
    TimeProviderBase& tp_;

    // order of RK solver
    const int ord_;
    // number of stages
    const int stages_;
    // init flag
    bool initialized_;
  };

  /** @} **/

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_EXPLICIT_HH
