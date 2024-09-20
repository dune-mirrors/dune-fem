#ifndef DUNE_FEM_SCHEMES_FEMSCHEME_HH
#define DUNE_FEM_SCHEMES_FEMSCHEME_HH

#include <dune/common/typeutilities.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/io/parameter/reader.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/solver/preconditionfunctionwrapper.hh>

namespace Dune
{
  namespace Fem
  {

// FemScheme
//----------

template < class Op, class DF, typename = void >
struct AddDirichletBC
{
  static const bool value = false;
  using DirichletBlockVector = void;
};
template < class Op, class DF>
struct AddDirichletBC<Op,DF,std::enable_if_t<std::is_void< decltype( std::declval<const Op>().
            setConstraints( std::declval<DF&>() ) )>::value > >
{
  static const bool value = true;
  using DirichletBlockVector = typename Op::DirichletBlockVector;
};

template< class Operator, class LinearInverseOperator,
          class InverseOperator = Dune::Fem::NewtonInverseOperator< typename Operator::JacobianOperatorType,
                                                                    LinearInverseOperator > >
class FemScheme
{
public:
  //! type of the mathematical model
  typedef typename Operator::ModelType ModelType;
  typedef typename Operator::DomainFunctionType DomainFunctionType;
  typedef typename Operator::RangeFunctionType  RangeFunctionType;
  typedef typename Operator::RangeFunctionType  DiscreteFunctionType;
  typedef Operator DifferentiableOperatorType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef LinearInverseOperator LinearInverseOperatorType;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;
  static_assert( std::is_same< typename DiscreteFunctionSpaceType::GridPartType, GridPartType >::value,
        "GridPart of Space has to be identical to GridPart of Model class" );

  //! type of underlying hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename Operator::JacobianOperatorType JacobianOperatorType;
  typedef typename Operator::JacobianOperatorType LinearOperatorType;

  // type of inverse operator (could be nonlinear or linear depending on the derived class)
  typedef InverseOperator InverseOperatorType;

  typedef typename InverseOperatorType::ErrorMeasureType ErrorMeasureType;

  typedef Dune::Fem::PreconditionerFunctionWrapper<
          typename LinearOperatorType::RangeFunctionType,
          typename LinearOperatorType::DomainFunctionType >  PreconditionerFunctionWrapperType;
  // std::function to represents the Python function passed as potential preconditioner
  typedef typename PreconditionerFunctionWrapperType::PreconditionerFunctionType  PreconditionerFunctionType ;

  typedef typename FunctionSpaceType::RangeType RangeType;
  static const int dimRange = FunctionSpaceType::dimRange;
  static constexpr bool addDirichletBC = AddDirichletBC<Operator,DomainFunctionType>::value;
  using DirichletBlockVector = typename AddDirichletBC<Operator,DomainFunctionType>::DirichletBlockVector;
  /*********************************************************/

  //! type of solver statistics reported (defined in inverseoperatorinterface.hh)
  typedef typename InverseOperatorType::SolverInfoType  SolverInfoType;

  //! constructor with one model
  FemScheme ( const DiscreteFunctionSpaceType &space, ModelType &model,
              const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : space_( space ),
    // the full discretized operator
    fullOpPtr_( new DifferentiableOperatorType(space, space, model, parameter) ),
    fullOperator_( *fullOpPtr_ ),
    // create inverse operator to invert the operator
    invOp_( parameter )
  {}

  //! constructor for derived classes (GalerkinScheme and MassLumpingScheme) with a list of models
  template < class... Models >
  FemScheme ( const DiscreteFunctionSpaceType &space,      // discrete function space
              const Dune::Fem::ParameterReader &parameter, // parameters
              Models &&... models )                        // list of models, could be more than one
  : space_( space ),
    // the full discretized operator
    fullOpPtr_( new DifferentiableOperatorType( space, space, std::forward< Models >( models )... )),
    fullOperator_( *fullOpPtr_ ),
    // create inverse operator to invert the operator
    invOp_( parameter )
  {}

  //! constructor for derived classes (LinearScheme and LinearizedScheme)
  FemScheme ( DifferentiableOperatorType& fullOp,
              const Dune::Fem::ParameterReader &parameter) // parameters
  : space_( fullOp.space() ),
    fullOpPtr_(), // empty here since fullOp reference is provided
    // the full discretized operator
    fullOperator_( fullOp ),
    // create inverse operator to invert the operator
    invOp_( parameter )
  {}

  const DifferentiableOperatorType &fullOperator() const { return fullOperator_; }
  DifferentiableOperatorType &fullOperator() { return fullOperator_; }

  template <typename O = DifferentiableOperatorType>
  auto setQuadratureOrders(unsigned int interior, unsigned int surface)
  -> Dune::void_t< decltype( std::declval< O >().setQuadratureOrders(0,0) ) >
  {
    fullOperator().setQuadratureOrders(interior,surface);
  }

  void setConstraints( DomainFunctionType &u ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().setConstraints( u );
  }
  void setConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().setConstraints( u,v );
  }
  template <class GridFunctionType>
  void setConstraints( const GridFunctionType &u, DiscreteFunctionType &v ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().setConstraints( u, v );
  }
  void setConstraints( const RangeType &value, DiscreteFunctionType &u ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().setConstraints( value, u );
  }
  void subConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().subConstraints( u, v );
  }
  void subConstraints( DiscreteFunctionType &v ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().subConstraints( v );
  }
  void addConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().addConstraints( u, v );
  }
  void addConstraints( DiscreteFunctionType &v ) const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      fullOperator().addConstraints( v );
  }
  const auto &dirichletBlocks() const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      return fullOperator().dirichletBlocks();
  }

  void operator() ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest ) const
  {
    fullOperator()( arg, dest );
  }
  template <class GridFunction>
  auto operator() ( const GridFunction &arg, DiscreteFunctionType &dest ) const
  -> Dune::void_t<decltype(std::declval<const Operator&>()(arg,dest))>
  {
    fullOperator()( arg, dest );
  }
  void setErrorMeasure(ErrorMeasureType &errorMeasure) const
  {
    invOp_.setErrorMeasure(errorMeasure);
  }

  SolverInfoType solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution) const
  {
    invOp_.bind(fullOperator());
    _solve(rhs,solution);
    invOp_.unbind();
    return invOp_.info();
  }
  SolverInfoType solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution, const PreconditionerFunctionType& p) const
  {
    PreconditionerFunctionWrapperType pre( p );
    invOp_.bind(fullOperator(), pre);
    _solve(rhs,solution);
    invOp_.unbind();
    return invOp_.info();
  }
  SolverInfoType solve ( DiscreteFunctionType &solution ) const
  {
    DiscreteFunctionType zero( solution );
    zero.clear();
    return solve(zero,solution);
  }
  SolverInfoType solve ( DiscreteFunctionType &solution, const PreconditionerFunctionType& p ) const
  {
    DiscreteFunctionType zero( solution );
    zero.clear();
    return solve(zero,solution,p);
  }

  template< class GridFunction, std::enable_if_t<
        std::is_same< decltype(
          std::declval< const DifferentiableOperatorType >().jacobian(
              std::declval< const GridFunction& >(), std::declval< JacobianOperatorType& >()
            )
          ), void >::value, int> i = 0
    >
  void jacobian( const GridFunction &ubar, JacobianOperatorType &linOp ) const
  {
    fullOperator().jacobian(ubar, linOp);
  }

  const GridPartType &gridPart () const { return space().gridPart(); }
  const DiscreteFunctionSpaceType &space( ) const { return space_; }

  const ModelType &model() const
  {
    return fullOperator().model();
  }
  ModelType &model()
  {
    return fullOperator().model();
  }
protected:
  SolverInfoType _solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution) const
  {                                     // setup for Newton scheme
    setConstraints(solution);           // sol=g     on bnd
    addConstraints(rhs,solution);       // sol=g+rhs on bnd
    invOp_( rhs, solution );
    return invOp_.info();
  }

  const DiscreteFunctionSpaceType &space_;   // discrete function space
  std::shared_ptr< DifferentiableOperatorType > fullOpPtr_; // full operator object storage
  DifferentiableOperatorType& fullOperator_; // reference to fullOperator (could be provided by derived class)
  mutable InverseOperatorType invOp_;        // non linear solver
};
} // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_SCHEMES_FEMSCHEME_HH
