#ifndef DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH
#define DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH

#include <cassert>
#include <iosfwd>
#include <string>

#include <dune/common/documentation.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/nullptr.hh>

#include <dune/fem/common/tupletypetraits.hh>
#include <dune/fem/common/tupleutility.hh>
#include <dune/fem/pass/common/localfunctiontuple.hh>
#include <dune/fem/pass/common/pass.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declarations
    // -----------------------------

    template< class DiscreteModel, class LocalOperator, class PreviousPass, int id >
    struct ApplyLocalOperatorPass;
    template< class DiscreteModel, class Argument, class PassIds, class Selector >
    class ApplyLocalOperatorDiscreteModelCaller;



    // ApplyLocalOperatorDiscreteModel
    // -------------------------------

    /** \brief Sample class layout of discrete models as expected by 
     *         ApplyLocalOperatorPass.
     *
     *  The discrete model determines the local function to be passed
     *  to the local operator in ApplyLocalOperatorPass.
     *
     */
    template< class DiscreteFunction >
    struct ApplyLocalOperatorDiscreteModel
    {
      typedef DiscreteFunction DiscreteFunctionType;
      //! \biref type of discrete function space type
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      //! \brief function space type
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief entity type
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
      //! \brief local coordinate type
      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

      //! \brief tuple of integral_constants containing pass ids
      typedef ImplementationDefined Selector;

      //! \brief set time
      void setTime ( double time );

      //! \brief return time
      double time () const;

      //! \brief return order of local function to be passed to local operator
      int order ( const EntityType &entity ) const;

      //! \brief evaluate method of local function 
      template< class ArgumentTuple >
      void evaluate ( const EntityType &entity,
                      const LocalCoordinateType &x,
                      const ArgumentTuple &tuple,
                      RangeType &value ) const;

      //! \brief jacobian method of local function
      template< class JacobianTuple >
      void jacobian ( const EntityType &entity,
                      const LocalCoordinateType &x,
                      const JacobianTuple &tuple,
                      JacobianRangeType &jacobian ) const;

      //! \brief hessian method of local function
      template< class HessianTuple >
      void hessian ( const EntityType &entity,
                     const LocalCoordinateType &x,
                     const HessianTuple &tuple,
                     HessianRangeType &hessian ) const;
    };



    // ApplyLocalOperatorPass
    // ----------------------

    /** \brief A pass implementation allowing for user defined local
     *         operation to be applied.
     *
     *  Local pass that evaluates a local function from given discrete
     *  model. The resulting local function is passed to the template 
     *  argument local operator. The local operator writes to the 
     *  destination of this pass.
     *
     *  The local operator must have the following form: 
     * \code
  struct LocalOperator
  {
    template< class LocalFunction, class LocalDofVector >
    void operator() ( const LocalFunction &localFunction, LocalDofVector &dofs ) const;
  };
     * \endcode
     *
     *  \tparam  DiscreteModel  discrete model
     *  \tparam  LocalOperator  local operator
     *  \tparam  PreviousPass   type of previous pass
     *  \tparam  id             pass id
     * 
     */
    template< class DiscreteModel, class LocalOperator, class PreviousPass, int id >
    class ApplyLocalOperatorPass
    : public Dune::Fem::LocalPass< DiscreteModel, PreviousPass, id >
    {
      typedef ApplyLocalOperatorPass< DiscreteModel, LocalOperator, PreviousPass, id > ThisType;
      typedef Dune::Fem::LocalPass< DiscreteModel, PreviousPass, id > BaseType;

    public:
      //! pass ids up to here (tuple of integral constants)
      typedef typename BaseType::PassIds PassIds;

      //! \brief type of discrete model
      typedef DiscreteModel DiscreteModelType;
      //! \brief type of local operator
      typedef LocalOperator LocalOperatorType;

      //! \brief argument type
      typedef typename BaseType::ArgumentType ArgumentType;
      //! \brief destination type
      typedef typename BaseType::DestinationType DestinationType;

      //! \brief discrete function space type
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      //! \brief entity type
      typedef typename BaseType::EntityType EntityType;

    private:
      typedef typename DiscreteModelType::Selector Selector;
      typedef ApplyLocalOperatorDiscreteModelCaller< DiscreteModelType, ArgumentType, PassIds, Selector > DiscreteModelCallerType;
      struct LocalFunction;

    public:
      using BaseType::space;
      using BaseType::time;

      ApplyLocalOperatorPass ( DiscreteModelType &discreteModel,
                               const LocalOperator &localOperator,
                               PreviousPass &pass,
                               const DiscreteFunctionSpaceType &space,
                               std::string passName = "ApplyLocalOperatorPass" )
      : BaseType( pass, space, passName ),
        discreteModel_( discreteModel ),
        localOperator_( localOperator ),
        arg_( nullptr ),
        dest_( nullptr ),
        caller_( nullptr ),
        localFunction_( nullptr )
      {}

      void printTexInfo (std::ostream &ostream ) const
      {
        BaseType::printTexInfo( ostream );
        ostream << "ApplyLocalOperatorPass: "
                << "\\\\ \n";
      }

    protected:
      void prepare ( const ArgumentType &arg, DestinationType &dest ) const
      {
        arg_ = const_cast< ArgumentType* >( &arg );
        dest_ = &dest;

        caller_ = new DiscreteModelCallerType( *arg_, discreteModel_ );
        assert( caller_ );
        caller_->setTime( time() );

        localFunction_ = new LocalFunction( *caller_ );
        assert( localFunction_ );
      }

      void finalize ( const ArgumentType &arg, DestinationType &dest ) const
      {
        space().communicate( dest );

        if( localFunction_ )
          delete localFunction_;

        if( caller_ )
          delete caller_;
        caller_ = nullptr;
      }

      void applyLocal ( const EntityType &entity ) const
      {
        localFunction_->init( entity );
        typename DestinationType::LocalFunctionType dofs = dest_->localFunction( entity );
        localOperator_( *localFunction_, dofs );
      }

    private:
      DiscreteModelType &discreteModel_;
      const LocalOperatorType &localOperator_;

      mutable ArgumentType *arg_;
      mutable DestinationType *dest_;
      mutable DiscreteModelCallerType *caller_;
      mutable LocalFunction *localFunction_;
    };



    // Implementation of ApplyLocalOperatorPass::LocalFunction
    // -------------------------------------------------------

    template< class DiscreteModel, class LocalOperator, class PreviousPass, int id >
    struct ApplyLocalOperatorPass< DiscreteModel, LocalOperator, PreviousPass, id >::LocalFunction
    {
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      static const int dimDomain = FunctionSpaceType::dimDomain;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef typename ApplyLocalOperatorPass< DiscreteModel, LocalOperator, PreviousPass, id >::EntityType EntityType;
      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

      LocalFunction ( DiscreteModelCallerType &caller )
      : caller_( caller )
      {}

      int order () const { return caller().order(); }

      const EntityType &entity () const { return caller_.entity(); }

      void init ( const EntityType &entity ) { caller_.setEntity( entity ); }

      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        caller().evaluate( x, value );
      }

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        caller().jacobian( x, jacobian );
      }

      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        caller().hessian( x, hessian );
      }

      template< class Quadrature, class RangeVector >
      void evaluateQuadrature( const Quadrature &quadrature, RangeVector &values ) const
      {
        caller().evaluateQuadrature( quadrature, values );
      }

    private:
      const DiscreteModelCallerType &caller () const { return caller_; }

      DiscreteModelCallerType &caller_;
    };



    // ApplyLocalOperatorDiscreteModelCaller
    // -------------------------------------

    template< class DiscreteModel, class Argument, class PassIds, class SelectorTuple >
    class ApplyLocalOperatorDiscreteModelCaller
    {
      typedef ApplyLocalOperatorDiscreteModelCaller< DiscreteModel, Argument, PassIds, SelectorTuple > ThisType;

    public:
      //! \brief discrete model type
      typedef DiscreteModel DiscreteModelType;
      //! \brief total argument type
      typedef Argument ArgumentType;
      //! \brief selector
      typedef SelectorTuple Selector;

      //! \brief entity type
      typedef typename DiscreteModelType::EntityType EntityType;

      //! \brief function space type
      typedef typename DiscreteModelType::FunctionSpaceType FunctionSpaceType;

      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    protected:
      typedef Dune::MakeSubTuple< ArgumentType, typename Dune::FirstTypeIndexTuple< PassIds, Selector >::type > FilterType;
      typedef typename FilterType::type DiscreteFunctionPointerTupleType;
      typedef typename TupleTypeTraits< DiscreteFunctionPointerTupleType >::PointeeTupleType DiscreteFunctionTupleType;

      typedef LocalFunctionTuple< DiscreteFunctionTupleType, EntityType > LocalFunctionTupleType;

      typedef typename LocalFunctionTupleType::RangeTupleType RangeTupleType;
      typedef typename LocalFunctionTupleType::JacobianRangeTupleType JacobianRangeTupleType;
      typedef typename LocalFunctionTupleType::HessianRangeTupleType HessianRangeTupleType;

    public:
      ApplyLocalOperatorDiscreteModelCaller ( ArgumentType &argument, DiscreteModel &discreteModel )
      : discreteModel_( discreteModel ),
        discreteFunctionPointerTuple_( FilterType::apply( argument ) ),
        localFunctionTuple_( DereferenceTuple< DiscreteFunctionPointerTupleType >::apply( discreteFunctionPointerTuple_ ) )
      {}
      
      //! \brief update local functions
      void setEntity ( const EntityType &entity )
      {
        localFunctionTuple().init( entity );
      }

      //! \brief set value to be returned by method time()
      void setTime ( double time ) // DUNE_DEPRECATED
      {
        discreteModel().setTime( time );
      }

      //!  \brief return time
      double time () const // DUNE_DEPRECATED
      {
        return discreteModel().time();
      }

      int order () const
      {
        return discreteModel().order( entity() );
      }

      //! \brief return entity
      const EntityType &entity () const { return localFunctionTuple().entity(); }

      //! \brief evalute local functions and pass values to discrete model
      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        localFunctionTuple().evaluate( x, values_ );
        discreteModel().evaluate( entity(), coordinate( x ), values_, value );
      }
      
      //! \brief evalute jacobians or local functions and pass values to discrete model
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        localFunctionTuple().evaluate( x, jacobians_ );
        discreteModel().jacobian( entity(), coordinate( x ), jacobians_, jacobian );
      }

      //! \brief evalute hessians of local functions and pass values to discrete model
      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        DUNE_THROW( Dune::NotImplemented, "Method hessian() not implemented yet" );
      }

      //! \brief calls method evaluate() for all quadrature points
      template< class Quadrature, class RangeVector >
      void evaluateQuadrature ( const Quadrature &quadrature, RangeVector &values ) const
      {
        assert( values.size() >= quadrature.nop() );
        const int nop = quadrature.nop();
        for( int qp = 0; qp < nop; ++qp )
          evaluate( entity(), quadrature[ qp ], values[ qp ] );
      }

    protected:
      DiscreteModelType &discreteModel () { return discreteModel_; }
      const DiscreteModelType &discreteModel () const { return discreteModel_; }

      LocalFunctionTupleType &localFunctionTuple () { return localFunctionTuple_; }
      const LocalFunctionTupleType &localFunctionTuple () const { return localFunctionTuple_; }

    private:
      ApplyLocalOperatorDiscreteModelCaller ( const ThisType & );
      ThisType &operator= ( const ThisType & );

      DiscreteModelType &discreteModel_;
      DiscreteFunctionPointerTupleType discreteFunctionPointerTuple_;
      LocalFunctionTupleType localFunctionTuple_;

      // caches for evaluating local functions
      mutable RangeTupleType values_;
      mutable JacobianRangeTupleType jacobians_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH
