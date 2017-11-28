#ifndef DUNE_FEM_SOLVER_ODESOLVERINTERFACE_HH
#define DUNE_FEM_SOLVER_ODESOLVERINTERFACE_HH

#include <iostream>

//- dune-common includes
#include <dune/common/exceptions.hh>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

namespace DuneODE
{
  using namespace Dune;
  using namespace Fem;

  /** \brief Interface class for ODE Solver. */
  template <class DestinationImp>
  class OdeSolverInterface
  {
  protected:
    //! constructor
    OdeSolverInterface () {}

    struct Monitor
    {
      double odeSolveTime_;
      double operatorTime_;
      double error_;

      std::size_t numberOfElements_;

      int newtonIterations_;
      int linearSolverIterations_;
      int maxNewtonIterations_;
      int maxLinearSolverIterations_;
      int spaceOperatorCalls_;

      Monitor() { reset(); }

      // reset all counters
      void reset()
      {
        odeSolveTime_ = 0;
        operatorTime_ = 0;
        error_ = 0;
        numberOfElements_ = 0;
        newtonIterations_ = 0;
        linearSolverIterations_ = 0;
        maxNewtonIterations_ = 0;
        maxLinearSolverIterations_ = 0;
        spaceOperatorCalls_ = 0;
      }
    };

  public:
    //! monitor type
    typedef Monitor MonitorType;

    //! type of destination
    typedef DestinationImp DestinationType;

    //! destructor
    virtual ~OdeSolverInterface () {}

    /** \brief initialize solver
        \param[in] arg argument to apply internal operator once for intial time step estimate
    */
    virtual void initialize(const DestinationType& arg) = 0;

    /** \brief solve \f$\partial_t u = L(u)\f$ where \f$L\f$ is the internal operator.
        \param[in] u unknown to solve for
    */
    virtual void solve(DestinationType& u)
    {
      MonitorType monitor;
      solve( u, monitor );
    }

    /** \brief solve \f$\partial_t u = L(u)\f$ where \f$L\f$ is the internal operator.
        \param[in] u unknown to solve for
        \param[in] monitor Monitor to get some inside information
    */
    virtual void solve ( DestinationType &u, MonitorType &monitor ) = 0;

    /** \brief print description of ODE solver to out stream */
    virtual void description(std::ostream&) const = 0;
  };


} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_ODESOLVERINTERFACE_HH
