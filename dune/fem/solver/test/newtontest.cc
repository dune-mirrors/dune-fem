// ****************************************
//
// Solve a non-linear equation 0 = F(x)
// with a Newton method
//
// ****************************************

// standard includes
#include <config.h>
#include <iostream>

// dune includes
#include <dune/common/fvector.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/solver/inverseoperatorinterface.hh>
#include <dune/common/fmatrix.hh>
#include <dune/fem/io/parameter.hh>

static const int systemSize = 5;

// Data structure for our unknown: fieldvector of
// dimension N (N = number of ODEs) with some additional methods.
// In scalar case (N=1) we have nothing more
// than a simple, tuned up double.
template <int N>
class DummyFunction : public Dune::FieldVector<double, N> {
  typedef DummyFunction< N > ThisType;
private:
  struct SpaceDummy {
    typedef double DomainFieldType;
    typedef double RangeFieldType;
    int size () const { return N; }
  };
  typedef Dune::FieldVector<double, N> BaseType;

public:
  typedef BaseType StorageType;
  typedef SpaceDummy DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  const DiscreteFunctionSpaceType space_;

  DummyFunction(std::string, const DiscreteFunctionSpaceType& space,
                std::initializer_list<RangeFieldType> const &u)
    : BaseType(u), space_( space )
  {}

  DummyFunction(std::string, std::initializer_list<RangeFieldType> const &u)
    : BaseType(u), space_()
  {}

  const DiscreteFunctionSpaceType& space () const { return space_; }

  DummyFunction() : space_( DiscreteFunctionSpaceType() ) {
    clear();
  }

  void clear() {
    BaseType::operator=(0.);
  }

  void assign(const DummyFunction& other) {
    BaseType::operator=(other);
  }

  void axpy(RangeFieldType l, const DummyFunction& other) {
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

  double normSquaredDofs() const
  {
    return scalarProductDofs( *this );
  }

};

template<int N>
class MyLinearOperator
  : public Dune::Fem::AssembledOperator<DummyFunction<N>, DummyFunction<N> >
{
  typedef typename DummyFunction<N>::StorageType VectorStorageType;
public:
  typedef Dune::FieldMatrix<typename DummyFunction<N>::RangeFieldType, N, N> MatrixType;
  typedef DummyFunction<N> DomainFunctionType;
  typedef DummyFunction<N> RangeFunctionType;


  template<class... Args>
  MyLinearOperator(const Args & ...)
  {}


  MatrixType& matrix()
  {
    return matrix_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

  /** \brief application operator
   *
   *  \param[in]   u  argument discrete function
   *  \param[out]  w  destination discrete function
   *
   *  \note This method has to be implemented by all derived classes.
   */
  void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
  {
    matrix_.mv(u, static_cast<VectorStorageType&>(w));
  }


  private:
  MatrixType matrix_;

};

template<int N> class MyLinearInverseOperator;

template< int N >
struct MyInverseOperatorTraits
{
  typedef DummyFunction<N>    DiscreteFunctionType;
  typedef Dune::Fem::Operator< DiscreteFunctionType, DiscreteFunctionType > OperatorType;
  typedef OperatorType  PreconditionerType;

  typedef OperatorType AssembledOperatorType;
  typedef DiscreteFunctionType SolverDiscreteFunctionType ;

  typedef MyLinearInverseOperator< N >  InverseOperatorType;
  typedef Dune::Fem::SolverParameter SolverParameterType;
};

template<int N>
class MyLinearInverseOperator
  : public Dune::Fem::LinearOperator<DummyFunction<N>, DummyFunction<N> >
  //: public Dune::Fem::InverseOperatorInterface< MyInverseOperatorTraits< N > >
{
  typedef Dune::Fem::InverseOperatorInterface< MyInverseOperatorTraits< N > > BaseType;
 public:
  //using BaseType::bind;
  typedef typename BaseType::DomainFunctionType DomainFunctionType;
  typedef typename BaseType::RangeFunctionType RangeFunctionType;
  typedef MyLinearOperator<N> LinearOperatorType;
  typedef Dune::Fem::SolverParameter SolverParameterType;

  MyLinearInverseOperator(const SolverParameterType& parameter = SolverParameterType() )
   : parameter_ ( parameter )
  { }

  void bind ( const LinearOperatorType& jOp ) { jOp_ = &jOp; }
  void unbind () { jOp_ = nullptr; }

  void operator()(const DomainFunctionType& u, RangeFunctionType& w) const
  {
    (*jOp_).matrix().solve(static_cast<typename DomainFunctionType::StorageType &> (w),static_cast<const typename DomainFunctionType::StorageType &> (u));
  }

  int apply(const DomainFunctionType& u, RangeFunctionType& w) const
  {
    (*jOp_).matrix().solve(static_cast<typename DomainFunctionType::StorageType &> (w),static_cast<const typename DomainFunctionType::StorageType &> (u));
    return 1;
  }
  int iterations() const { return 1; }
  void setMaxIterations ( unsigned int ) {}

  SolverParameterType &parameter() const {
    return parameter_;
  }

 private:
  const LinearOperatorType * jOp_ = nullptr;
  mutable SolverParameterType parameter_;
};




template<int N>
class MyOperator
  : public Dune::Fem::DifferentiableOperator<MyLinearOperator<N> >
{
  typedef Dune::Fem::DifferentiableOperator<MyLinearOperator<N> > BaseType;
  const bool lineSearch_;
public:
  MyOperator( const bool lineSearch ) : lineSearch_( lineSearch ) {}

  typedef typename BaseType::DomainFunctionType DomainFunctionType;
  typedef typename BaseType::RangeFunctionType RangeFunctionType;
  typedef typename BaseType::JacobianOperatorType JacobianOperatorType;

  void operator()(const DomainFunctionType& u, RangeFunctionType& w) const
  {
    for (int i = 0; i < N; ++i)
    {
      if( lineSearch_ )
        w[i] = std::atan(u[i]);
      else
        w[i] = u[i]*u[i] - i;
    }
  }

  void jacobian(const DomainFunctionType& u, JacobianOperatorType& jOp) const
  {
    auto& matrix = jOp.matrix();

    matrix = 0.;
    for (int i = 0; i < N; ++i)
    {
      if( lineSearch_ )
        matrix[i][i] = 1./(1.+u[i]*u[i]);
      else
        matrix[i][i] = 2*u[i];
    }
  }

};

template <class FunctionType>
void run(FunctionType& sol, FunctionType& rhs, const bool lineSearch )
{
  typedef MyLinearOperator<systemSize> LinearOperatorType;
  typedef MyOperator<systemSize> OperatorType;
  typedef MyLinearInverseOperator<systemSize> LinearInverseOperatorType;
  typedef Dune::Fem::NewtonInverseOperator<LinearOperatorType, LinearInverseOperatorType> NewtonInverseOperatorType;

  OperatorType op( lineSearch );
  NewtonInverseOperatorType opInv;
  opInv.bind( op );

  //solve op(sol) = rhs for sol
  //with a Newton method
  //and initial value sol
  opInv(rhs, sol);

  std::cout << "Solution: " << sol << std::endl;
  std::cout << "Right Hand Side: " << rhs << std::endl;

  std::cout << "#Iterations: " << opInv.iterations() << std::endl;

  std::cout << "Converged: " << (opInv.converged() ? "yes" : "no") << std::endl;

  opInv.unbind();
}

int main( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  Dune::Fem::Parameter::append( argc, argv );
  if( argc == 1 )
    Dune::Fem::Parameter::append("parameter");

  const bool verbose = Dune::Fem::Parameter::verbose();
  (void) verbose;

  typedef DummyFunction<systemSize> FunctionType;

  // standard Newton
  {
    FunctionType sol("sol", { 1, 2, 3, 4, 7 }), rhs("rhs", {0,0,0,0,0});
    run( sol, rhs, false );
  }

  // Newton with line search
  {
    Dune::Fem::Parameter::append("fem.solver.newton.lineSearch","simple");
    FunctionType sol("sol", { 2, 2, 2, 2, 2 }), rhs("rhs", { -1.57, -1.5, 0, 1.5, 1.57 });
    run( sol, rhs, true );
  }

  return 0;
}
