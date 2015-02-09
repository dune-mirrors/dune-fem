// ****************************************
//
// Solve the system of ODEs dy/dt = F(y,t)
//                         F(y,0) = 0
//
// ****************************************
#undef ENABLE_MPI

// standard includes
#include <config.h>
#include <iostream>

// dune includes
#include <dune/common/fvector.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#include <dune/fem/solver/odesolver.hh>

#include <dune/fem/solver/rungekutta/explicit.hh>

#include <dune/fem/io/parameter.hh>
//#include <dune/fem/solver/rungekutta/implicit.hh>

static const int systemSize = 4;

// Data structure for our unknown: fieldvector of
// dimension N (N = number of ODEs) with some additional methods.
// In scalar case (N=1) we have nothing more
// than a simple, tuned up double.
template <int N>
class myDest : public Dune::FieldVector<double, N> {
  typedef myDest< N > ThisType;
private:
  struct SpaceDummy {
    int size () const { return N; }
  };
  typedef Dune::FieldVector<double, N> BaseType;

public:
  typedef double DomainFieldType;
  typedef double RangeFieldType;
  typedef SpaceDummy DiscreteFunctionSpaceType;

  myDest(std::string, const SpaceDummy&, const double* u = 0) {
    clear();
  }

  myDest() {
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

  double* leakPointer() { return &this->operator[](0); }
  const double* leakPointer() const { return &this->operator[](0); }

  double scalarProductDofs ( const ThisType &other ) const
  {
    double scp = 0;
    for( std::size_t i=0; i < N; ++i )
      scp += (*this)[ i ] * other[ i ];

    return scp;
  }

  double timeStepEstimate() const { return 0.01; }
};


// implement right hand side F(y,t)
// here: system of three ODEs
class myRHS : public Dune::Fem::SpaceOperatorInterface< myDest<systemSize> > {
public:
  myRHS() {
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
  }

  static DestinationType exact (const double t) {
    DestinationType exact;
    exact[ 0 ] = t * t;
    exact[ 1 ] = t * t * t;
    exact[ 2 ] = 0 ;
    exact[ 4 ] = -3.0/5.0 * std::exp( 5.0 * t ) + 3.0/5.0;
    return exact;
  }

  void setTime(const double time) {
    t_=time;
  }

private:
  SpaceType space_;
  double t_;
};


template <class OdeSolverType>
void solve(const bool verbose)
{
  typedef myRHS SpaceOperatorType;
  typedef SpaceOperatorType::DestinationType DestinationType;

  // problem data
  const double startTime = 0.0;
  const double endTime = 2.0;

  // options
  const double stepSize = 0.00001;
  const double cfl = 1.;
  const int order = 3;

  // create solver
  Dune::Fem::DefaultTimeProvider tp( startTime, cfl );
  SpaceOperatorType spaceOperator;
  OdeSolverType odeSolver( spaceOperator, tp, order );

  // initialize solution vector, same initial data for all components
  DestinationType U;
  U.clear(); // set initial data U(0) = 0

  // initialize odesolver
  odeSolver.initialize( U );

  // time loop
  for( tp.init(stepSize); tp.time() < endTime; tp.next(stepSize) ) {
    // do calculation
    odeSolver.solve(U);

    if( verbose )
    {
      // print out solution
      std::cout << "time = " << tp.time() << ", U = "
                << U
                << std::endl;
    }
  }

  // print out solution
  std::cout << "Result (t = " << tp.time() << ") U = " << U << std::endl;
  std::cout << "                 ex U = " << spaceOperator.exact( tp.time() ) << std::endl;
}


int main( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  Dune::Fem::Parameter::append( argc, argv );
  if( argc == 1 )
    Dune::Fem::Parameter::append("parameter");

  const bool verbose = Dune::Fem::Parameter::verbose();

  // types
  typedef myRHS SpaceOperatorType;
  typedef SpaceOperatorType::DestinationType DestinationType;

  /*
  {
    typedef DuneODE::ExplicitOdeSolver<DestinationType> OdeSolverType;
    solve< OdeSolverType > ( verbose );
  }

  */
  {
    typedef DuneODE::ImplicitOdeSolver<DestinationType> OdeSolverType;
    solve< OdeSolverType > ( verbose );
  }

  return 0;
}
