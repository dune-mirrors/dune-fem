// ****************************************
//
// Solve the system of ODEs dy/dt = F(y,t)
//                         F(y,0) = 0
//
// ****************************************

// standard includes
#include <config.h>
#include <iostream>

// dune includes
#include <dune/common/fvector.hh>
#include <dune/common/parallel/communication.hh>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#include <dune/fem/solver/odesolver.hh>

#include <dune/fem/solver/rungekutta/timestepcontrol.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/operator/dghelmholtz.hh>

#include <dune/fem/solver/rungekutta/explicit.hh>
#include <dune/fem/solver/rungekutta/implicit.hh>
#include <dune/fem/solver/rungekutta/row.hh>

#include <dune/fem/io/parameter.hh>

static const int systemSize = 5;

// Data structure for our unknown: fieldvector of
// dimension N (N = number of ODEs) with some additional methods.
// In scalar case (N=1) we have nothing more
// than a simple, tuned up double.
template <int N>
class myDest : public Dune::FieldVector<double, N> {
  typedef myDest< N > ThisType;
private:
  struct GridPartDummy
  {
    typedef Dune::Communication< Dune::No_Comm > CommunicationType;

    const CommunicationType &comm () const { return comm_; }

    CommunicationType comm_;
  };

  struct SpaceDummy {
    typedef GridPartDummy GridPartType;
    typedef std::array<int, 1 > AuxiliaryDofsType;

    SpaceDummy () : auxiliaryDofs_() { auxiliaryDofs_[ 0 ] = N; }

    int size () const { return N; }

    const GridPartType &gridPart () const { return gridPart_; }
    const AuxiliaryDofsType &auxiliaryDofs () const { return auxiliaryDofs_; }

    AuxiliaryDofsType auxiliaryDofs_;
    GridPartType gridPart_;
  };
  typedef Dune::FieldVector<double, N> BaseType;

public:
  typedef double DomainFieldType;
  typedef double RangeFieldType;
  typedef SpaceDummy DiscreteFunctionSpaceType;
  const DiscreteFunctionSpaceType space_;

  typedef typename BaseType::const_iterator ConstDofIteratorType;
  typedef typename BaseType::iterator DofIteratorType;

  myDest(std::string, const DiscreteFunctionSpaceType& space, const double* u = 0)
   : space_( space )
  {
    if( u )
    {
      for( int i=0; i<N; ++i )
      {
        (*this)[i] = u[i];
      }
    }
    else
      clear();
  }

  const DiscreteFunctionSpaceType& space () const { return space_; }

  myDest() : space_( DiscreteFunctionSpaceType() ) {
    clear();
  }

  void clear() {
    BaseType::operator=(0.);
  }

  void assign(const myDest& other) {
    BaseType::operator=(other);
  }

  void axpy(RangeFieldType l, const myDest& other) {
    BaseType::axpy(l, other);
  }

  double operator()(int i) const {
    if (i<0 || i>N-1) {
      std::cout << "ERROR: Accessing element " << i << std::endl;
    }
    return (*this)[i];
  }

  myDest& dofVector() { return *this; }
  const myDest& dofVector() const { return *this; }

  ConstDofIteratorType dbegin () const { return static_cast< const BaseType & >( *this ).begin(); }
  DofIteratorType dbegin () { return static_cast< BaseType & >( *this ).begin(); }
  ConstDofIteratorType dend () const { return static_cast< const BaseType & >( *this ).end(); }
  DofIteratorType dend () { return static_cast< BaseType & >( *this ).end(); }

  double* leakPointer() { return &this->operator[](0); }
  const double* leakPointer() const { return &this->operator[](0); }

  double scalarProductDofs ( const ThisType &other ) const
  {
    double scp = 0;
    for( std::size_t i=0; i < N; ++i )
      scp += (*this)[ i ] * other[ i ];

    return scp;
  }

  double normSquaredDofs() const
  {
    return scalarProductDofs( *this );
  }

};

// implement right hand side F(y,t)
// here: system of three ODEs
class myRHS : public Dune::Fem::SpaceOperatorInterface< myDest<systemSize> > {
public:
  typedef myDest<systemSize> DestinationType ;

  typedef myRHS PreconditionOperatorType;

  myRHS() : t_(0.0) {
  }

  const SpaceType& space() const {
    return space_;
  }

  void operator()(const DestinationType& X,
                  DestinationType& Y) const
  {
    const double* x = X.leakPointer();
    double* y = Y.leakPointer();
    this->operator()( x, y );
  }

  void operator()(const double* x,
                  double* y) const
  {
    y[0] = 2.0*t_;
    y[1] = 3.0*t_*t_;
    y[2] = x[2];
    y[3] = 5.0 * x[ 3 ] - 3.0;
    y[4] = std::cos( t_ );
  }

  static DestinationType exact (const double t) {
    DestinationType exact;
    exact[ 0 ] = t * t;
    exact[ 1 ] = t * t * t;
    exact[ 2 ] = 0 ;
    exact[ 3 ] = -3.0/5.0 * std::exp( 5.0 * t ) + 3.0/5.0;
    exact[ 4 ] = std::sin( t );
    return exact;
  }

  void setTime(const double time) {
    t_=time;
  }

  double timeStepEstimate() const { return 0.01; }

  PreconditionOperatorType* preconditioner() const
  {
    return 0;
  }

private:
  SpaceType space_;
  double t_;
};


#if 0
namespace Dune
{
  namespace Fem
  {

    namespace LinearSolver
    {
      template <>
      class OperatorAdapter< myDest<systemSize>, myDest<systemSize> >
      : public LinearSolver::FunctionIF<double>
      {
      public:
        typedef myDest<systemSize> DomainFunctionType;
        typedef DomainFunctionType RangeFunctionType;
        typedef OperatorAdapter< DomainFunctionType, RangeFunctionType> ThisType;

      public:
        typedef Operator< DomainFunctionType, RangeFunctionType > OperatorType;

        typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

        OperatorAdapter ( const OperatorType &op, const DomainFunctionSpaceType &domainSpace, const RangeFunctionSpaceType &rangeSpace )
        : operator_( op ),
          domainSpace_( domainSpace ),
          rangeSpace_( rangeSpace )
        {}

        void operator() ( const double *u, double *w)
        {
          DomainFunctionType uFunction( "OperatorAdapter::u", domainSpace_, u );
          RangeFunctionType  wFunction( "OperatorAdapter::w", rangeSpace_, w );
          operator_( uFunction, wFunction );
          // copy result back
          const int dim = size();
          for( int i=0; i<dim; ++ i )
            w[ i ] = wFunction[ i ];
        }

        int size() const
        {
          return domainSpace_.size();
        }

      private:
        const OperatorType &operator_;
        const DomainFunctionSpaceType &domainSpace_;
        const RangeFunctionSpaceType &rangeSpace_;
      };
    }
  }
}
#endif

template <class OdeFactory>
void solve(OdeFactory factory, const bool verbose)
{
  typedef myRHS SpaceOperatorType;
  typedef SpaceOperatorType::DestinationType DestinationType;
  SpaceOperatorType spaceOperator;

  // problem data
  const double startTime = 0.0;
  const double endTime = 2.0;

  // options
  const double stepSize = Dune::Fem::Parameter::getValue("fixedTimeStep", double(spaceOperator.timeStepEstimate()/2.)) ;
  std::cout << "dt = " << stepSize << std::endl;
  const double cfl = 1.;
  const int order = Dune::Fem::Parameter::getValue("fem.ode.order", int(2));

  // create solver
  Dune::Fem::TimeProvider<> tp( startTime, cfl );
  typedef typename OdeFactory :: OdeSolverType OdeSolverType;
  std::unique_ptr< OdeSolverType > odeSolver;
  odeSolver.reset( factory.create( spaceOperator, tp, order ) );

  // initialize solution vector, same initial data for all components
  DestinationType U;
  U.clear(); // set initial data U(0) = 0

  // initialize odesolver
  odeSolver->initialize( U );

  // time loop
  for( tp.init(stepSize); tp.time() < endTime; tp.next(stepSize) ) {
    // do calculation
    odeSolver->solve(U);

    if( verbose )
    {
      // print out solution
      std::cout << "time = " << tp.time() << ", U = "
                << U
                << std::endl;
    }
  }

  // print out solution
  std::cout << "Result (t = " << std::setw(6) << tp.time() << ") U = " << U << std::endl;
  DestinationType exact = spaceOperator.exact( tp.time() );
  std::cout << "              exact U = " << exact << std::endl;
  exact -= U;
  exact /= U.two_norm();

  auto twonorm = exact.two_norm();
  std::cout << "Two norm: " << twonorm << std::endl;
  if( std::to_string( twonorm ) == "nan" || twonorm > 1e-2 )
  {
    DUNE_THROW(Dune::InvalidStateException,"ERROR: ode solver did not converge!");
  }
}

template <class OdeSolver>
struct SimpleFactory
{
  typedef OdeSolver  OdeSolverType;
  template <class SpaceOperatorType, class TimeProvider>
  OdeSolverType* create( SpaceOperatorType& op, TimeProvider& tp, const int order )
  {
    return new OdeSolverType( op, tp, order );
  }
};

template <class SpaceOperator, template <class,class,class> class Solver >
struct ImplicitRKFactory
{
  typedef SpaceOperator  SpaceOperatorType;
  typedef typename SpaceOperatorType::DestinationType DestinationType;

  typedef Dune::Fem::KrylovInverseOperator< DestinationType >               LinearInverseOperatorType;
  typedef DuneODE::ImplicitRungeKuttaTimeStepControl                        TimeStepControlType;

  typedef Dune::Fem::DGHelmholtzOperator< SpaceOperatorType >               HelmholtzOperatorType;

  typedef Dune::Fem::NewtonInverseOperator< typename HelmholtzOperatorType::JacobianOperatorType,
                                            LinearInverseOperatorType >     NonlinearInverseOperatorType;

  // either ImplicitRungeKuttaSolver or ROWRungeKuttaSolver
  typedef Solver< HelmholtzOperatorType, NonlinearInverseOperatorType, TimeStepControlType >   OdeSolverType;

  std::unique_ptr< HelmholtzOperatorType > helmOp_;

  template <class TimeProvider>
  OdeSolverType* create( SpaceOperatorType& op, TimeProvider& tp, int order )
  {
    helmOp_.reset( new HelmholtzOperatorType( op ) );
    return new OdeSolverType( *helmOp_, tp, order );
  }
};

int main( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  Dune::Fem::Parameter::append( argc, argv );
  if( argc == 1 )
    Dune::Fem::Parameter::append("parameter");

  const bool verbose = Dune::Fem::Parameter::verbose();

  // types
  typedef myRHS SpaceOperatorType;
  typedef SpaceOperatorType::DestinationType DestinationType;

  // explicit RungeKutta (dune expl)
  {
    std::cout << "Dune-fem explicit rungekutta" << std::endl;
    typedef DuneODE::ExplicitRungeKuttaSolver<DestinationType> OdeSolverType;
    solve( SimpleFactory< OdeSolverType >(), verbose );
  }

  // implicit RungeKutta (dune impl)
  {
    std::cout << "Dune-fem implicit rungekutta" << std::endl;
    solve( ImplicitRKFactory< SpaceOperatorType, DuneODE::ImplicitRungeKuttaSolver > (), verbose );
  }

  // row RungeKutta (dune row)
  {
    std::cout << "Dune-fem row rungekutta" << std::endl;
    solve( ImplicitRKFactory< SpaceOperatorType, DuneODE::ROWRungeKuttaSolver > (), verbose );
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;

}
