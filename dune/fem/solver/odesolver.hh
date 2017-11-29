#ifndef RUNGEKUTTA_ODE_SOLVER_HH
#define RUNGEKUTTA_ODE_SOLVER_HH

//- system includes
#include <iostream>
#include <cmath>
#include <vector>
#include <pthread.h>
#include <cassert>
#include <sys/times.h>

//- dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/operator/dghelmholtz.hh>

#include <dune/fem/solver/odesolverinterface.hh>
#include <dune/fem/solver/timeprovider.hh>

#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>

#include <dune/fem/solver/rungekutta/timestepcontrol.hh>
#include <dune/fem/solver/rungekutta/explicit.hh>
#include <dune/fem/solver/rungekutta/implicit.hh>
#include <dune/fem/solver/rungekutta/semiimplicit.hh>

namespace DuneODE
{
  using namespace Dune;
  using namespace Fem;
  using namespace std;

  using ODEParameters = ImplicitRungeKuttaSolverParameters ;

  /**
     @ingroup ODESolver
     @{
   **/

  ///////////////////////////////////////////////////////
  //
  //  --ExplicitOdeSolver
  //
  ///////////////////////////////////////////////////////
  template< class Destination >
  using ExplicitOdeSolver = ExplicitRungeKuttaSolver< Destination >;

  ///////////////////////////////////////////////////////
  //
  //  --ImplicitOdeSolver
  //
  ///////////////////////////////////////////////////////
  template <class Destination>
  using ParDGHelmholtz = DGHelmholtzOperator< SpaceOperatorInterface<Destination> >;

  template <class Destination>
  using ParDGNewtonInverse = NewtonInverseOperator<
          typename ParDGHelmholtz< Destination > :: JacobianOperatorType,
          KrylovInverseOperator< Destination > >;

  template<class Destination>
  class ImplicitOdeSolver
    : public ImplicitRungeKuttaSolver<
          ParDGHelmholtz< Destination >, ParDGNewtonInverse< Destination > >
  {
  public:
    typedef SpaceOperatorInterface<Destination> OperatorType;
    typedef ParDGHelmholtz< Destination >       HelmholtzOperatorType;
    typedef ImplicitRungeKuttaSolver< HelmholtzOperatorType, ParDGNewtonInverse< Destination > > BaseType;

  protected:
    using BaseType :: helmholtzOp_;
    std::unique_ptr< HelmholtzOperatorType > helmOpPtr_;

    HelmholtzOperatorType* createHelmholtzOperator( OperatorType& op )
    {
      return new HelmholtzOperatorType( op ) ;
    }

  public:
    ImplicitOdeSolver( OperatorType& op, TimeProviderBase& tp, int order,
                       const ParameterReader &parameter = Parameter::container() )
      : BaseType( *createHelmholtzOperator( op ), tp, order, parameter )
    {
      // store pointer for later removal
      helmOpPtr_.reset( &helmholtzOp_ );
    }
  };

  ///////////////////////////////////////////////////////
  //
  //  --SemiImplicitOdeSolver
  //
  ///////////////////////////////////////////////////////
  template<class Destination>
  class SemiImplicitOdeSolver
    : public SemiImplicitRungeKuttaSolver<
          SpaceOperatorInterface<Destination>,
          ParDGHelmholtz< Destination >, ParDGNewtonInverse< Destination > >
  {
  public:
    typedef SpaceOperatorInterface<Destination> OperatorType;
    typedef ParDGHelmholtz< Destination >       HelmholtzOperatorType;
    typedef SemiImplicitRungeKuttaSolver< OperatorType, HelmholtzOperatorType, ParDGNewtonInverse< Destination > > BaseType;

  protected:
    using BaseType :: helmholtzOp_;
    std::unique_ptr< HelmholtzOperatorType > helmOpPtr_;

    HelmholtzOperatorType* createHelmholtzOperator( OperatorType& op )
    {
      return new HelmholtzOperatorType( op ) ;
    }

  public:
    SemiImplicitOdeSolver( OperatorType& explOp, OperatorType& implOp, TimeProviderBase& tp, int order,
                       const ParameterReader &parameter = Parameter::container() )
      : BaseType( explOp, *createHelmholtzOperator( implOp ), tp, order, parameter )
    {
      // store pointer for later removal
      helmOpPtr_.reset( &helmholtzOp_ );
    }
  };

  /**
   @}
  **/

} // namespace DuneODE

#endif // #ifndef RUNGEKUTTA_ODE_SOLVER_HH
