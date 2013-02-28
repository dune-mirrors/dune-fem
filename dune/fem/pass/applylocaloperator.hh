#ifndef DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH
#define DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH

#include <cassert>
#include <iosfwd>
#include <string>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/nullptr.hh>

#include <dune/fem/pass/common/filter.hh>
#include <dune/fem/pass/common/localfunctiontuple.hh>
#include <dune/fem/pass/common/pointertuple.hh>
#include <dune/fem/pass/common/selector.hh>
#include <dune/fem/pass/common/tupletypetraits.hh>
#include <dune/fem/pass/pass.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declarations
    // -----------------------------

    template< class DiscreteModel, class LocalOperator, class PreviousPass, int id >
    struct ApplyLocalOperatorPass;
    template< class DiscreteModel, class Argument, class Pass, class Selector >
    class ApplyLocalOperatorDiscreteModelCaller;



    // ApplyLocalOperatorDiscreteModel
    // -------------------------------

    template< class Traits,
              int N1 = -1,
              int N2 = -1,
              int N3 = -1,
              int N4 = -1,
              int N5 = -1,
              int N6 = -1,
              int N7 = -1,
              int N8 = -1,
              int N9 = -1
            >
    struct ApplyLocalOperatorDiscreteModel
    {
      typedef typename Traits::DiscreteModelType DiscreteModelType;
      typedef typename Traits::DestinationType DiscreteFunctionType;
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

      typedef typename Dune::Fem::Selector< N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9 >::Type Selector;

      void setTime ( double time )
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().setTime( time ) );
        return asImp().setTime( time );
      }

      double time () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().time() );
        return asImp().time();
      }

      int order ( const EntityType &entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().order( entity ) );
        return asImp().order( entity );
      }

      template< class ArgumentTuple >
      void evaluate ( const EntityType &entity,
                      const LocalCoordinateType &x,
                      const ArgumentTuple &tuple,
                      RangeType &value ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().evaluate( entity, x, tuple, value ) );
        return asImp().evaluate( entity, x, tuple, value );
      }

      template< class JacobianTuple >
      void jacobian ( const EntityType &entity,
                      const LocalCoordinateType &x,
                      const JacobianTuple &tuple,
                      JacobianRangeType &jacobian ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().jacobian( entity, x, tuple, jacobian ) );
        return asImp().evaluate( entity, x, tuple, jacobian );
      }

      template< class HessianTuple >
      void hessian ( const EntityType &entity,
                     const LocalCoordinateType &x,
                     const HessianTuple &tuple,
                     HessianRangeType &hessian ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().evaluate( entity, x, tuple, hessian ) );
        return asImp().evaluate( entity, x, tuple, hessian );
      }

    protected:
      DiscreteModelType &asImp () { return static_cast< DiscreteModelType & >( *this ); }
      const DiscreteModelType &asImp () const { return static_cast< const DiscreteModelType & >( *this ); }
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
    struct ApplyLocalOperatorPass 
    : public Dune::Fem::LocalPass< DiscreteModel, PreviousPass, id >
    {
      typedef ApplyLocalOperatorPass< DiscreteModel, LocalOperator, PreviousPass, id > ThisType;
      typedef Dune::Fem::LocalPass< DiscreteModel, PreviousPass, id > BaseType;

    public:
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
      typedef ApplyLocalOperatorDiscreteModelCaller< DiscreteModelType, ArgumentType, ThisType, Selector > DiscreteModelCallerType;
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
        caller_( discreteModel ),
        localFunction_( caller_ ),
        localOperator_( localOperator ),
        arg_( nullptr ),
        dest_( nullptr )
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
        caller_.prepare( *arg_ );
        caller_.setTime( time() );
      }

      void finalize ( const ArgumentType &arg, DestinationType &dest ) const
      {
        space().communicate( dest );
        caller_.finalize();
      }

      void applyLocal ( const EntityType &entity ) const
      {
        // initialize argument for local operator
        LocalFunction &arg = localFunction_;
        arg.init( entity );
        
        // get destination for local operator
        typename DestinationType::LocalFunctionType dest = dest_->localFunction( entity );

        // apply local operator
        localOperator_( arg, dest );
      }

    private:
      mutable DiscreteModelCallerType caller_;
      mutable LocalFunction localFunction_;
      const LocalOperatorType &localOperator_;

      mutable ArgumentType *arg_;
      mutable DestinationType *dest_;
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
      : caller_( caller ),
        entity_( nullptr )
      {}

      int order () const { return caller().order( entity() ); }

      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      void init ( const EntityType &entity )
      {
        entity_ = &entity;
        caller_.setEntity( entity );
      }

      template< class PointType >
      void evaluate ( const PointType &x, RangeType &value ) const
      {
        caller().evaluate( entity(), x, value );
      }

      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &jacobian ) const
      {
        caller().jacobian( entity(), x, jacobian );
      }

      template< class PointType >
      void hessian ( const PointType &x, HessianRangeType &hessian ) const
      {
        caller().hessian( entity(), x, hessian );
      }

    private:
      const DiscreteModelCallerType &caller () const { return caller_; }

      DiscreteModelCallerType &caller_;
      const EntityType *entity_;
    };



    // ApplyLocalOperatorDiscreteModelCaller
    // -------------------------------------

    template< class DiscreteModel, class Argument, class Pass, class SelectorTuple >
    class ApplyLocalOperatorDiscreteModelCaller
    {
      typedef ApplyLocalOperatorDiscreteModelCaller< DiscreteModel, Argument, Pass, SelectorTuple > ThisType;

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
      typedef Filter< ArgumentType, Pass, Selector > FilterType;
      typedef PointerTuple< typename FilterType::ResultType > DiscreteFunctionPointerTupleType;
      typedef typename DiscreteFunctionPointerTupleType::ElementType DiscreteFunctionTupleType;

      typedef LocalFunctionTuple< DiscreteFunctionTupleType, EntityType > LocalFunctionTupleType;

      typedef typename LocalFunctionTupleType::RangeTupleType RangeTupleType;
      typedef typename LocalFunctionTupleType::JacobianRangeTupleType JacobianRangeTupleType;
      typedef typename LocalFunctionTupleType::HessianRangeTupleType HessianRangeTupleType;

    public:
      ApplyLocalOperatorDiscreteModelCaller ( DiscreteModel &discreteModel )
      : discreteModel_( discreteModel ),
        localFunctionTuple_( nullptr )
      {
        // we are going to use placement new
        buffer_ = new char[ sizeof( LocalFunctionTupleType ) ];
      }

      ~ApplyLocalOperatorDiscreteModelCaller ()
      {
        if( localFunctionTuple_ )
          localFunctionTuple_->~LocalFunctionTupleType();
        localFunctionTuple_ = nullptr;

        delete[] buffer_;
      }

      //! \brief prepare
      void prepare ( ArgumentType &argument )
      {
        discreteFunctionPointerTuple_ = DiscreteFunctionPointerTupleType( FilterType::apply( argument ) );
        localFunctionTuple_ = new( buffer_) LocalFunctionTupleType( *discreteFunctionPointerTuple_ );
      }

      //! \brief finalize
      void finalize ()
      {
        if( localFunctionTuple_ )
          localFunctionTuple_->~LocalFunctionTupleType();
        localFunctionTuple_ = nullptr;

        discreteFunctionPointerTuple_ = DiscreteFunctionPointerTupleType();
      }
      
      //! \brief update local functions
      void setEntity ( const EntityType &entity )
      {
        localFunctionTuple().setEntity( entity );
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

      int order ( const EntityType &entity ) const
      {
        return discreteModel().order( entity );
      }

      //! \brief evalute local functions and pass values to discrete model
      template< class Point >
      void evaluate ( const EntityType &entity,
                      const Point &x,
                      RangeType &value ) const
      {
        assert( &entity == &(localFunctionTuple().entity()) );
        localFunctionTuple().evaluate( x, values_ );
        discreteModel().evaluate( entity, coordinate( x ), values_, value );
      }
      
      //! \brief evalute jacobians or local functions and pass values to discrete model
      template< class Point >
      void jacobian ( const EntityType &entity,
                      const Point &x,
                      JacobianRangeType &jacobian ) const
      {
        assert( &entity == &(localFunctionTuple().entity()) );
        localFunctionTuple().evaluate( x, jacobians_ );
        discreteModel().jacobian( entity, coordinate( x ), jacobians_, jacobian );
      }

      //! \brief evalute hessians of local functions and pass values to discrete model
      template< class Point >
      void hessian ( const EntityType &entity,
                     const Point &x,
                     HessianRangeType &hessian ) const
      {
        DUNE_THROW( Dune::NotImplemented, "Method hessian() not implemented yet" );
      }

      //! \brief calls method evaluate() for all quadrature points
      template< class QuadratureType, class RangeVectorType >
      void evaluateQuadrature ( const EntityType &entity,
                                const QuadratureType &quadrature,
                                RangeVectorType &values ) const
      {
        const int nop = quadrature.nop();
        for( int qp = 0; qp < nop; ++qp )
          evaluate( entity, quadrature[ qp ], values[ qp ] );
      }

    protected:
      DiscreteModelType &discreteModel () { return discreteModel_; }

      const DiscreteModelType &discreteModel () const { return discreteModel_; }

      LocalFunctionTupleType &localFunctionTuple ()
      {
        assert( localFunctionTuple_ );
        return *localFunctionTuple_;
      }

      const LocalFunctionTupleType &localFunctionTuple () const
      {
        assert( localFunctionTuple_ );
        return *localFunctionTuple_;
      }

      static DiscreteFunctionPointerTupleType nullptrTuple ()
      {
        return Dune::NullPointerInitialiser< DiscreteFunctionTupleType >::apply();
      }

    private:
      ApplyLocalOperatorDiscreteModelCaller ( const ThisType & );
      ThisType &operator= ( const ThisType & );

      DiscreteModelType &discreteModel_;
      
      DiscreteFunctionPointerTupleType discreteFunctionPointerTuple_;

      LocalFunctionTupleType *localFunctionTuple_;
      char *buffer_;

      // caches for evaluating local functions
      mutable RangeTupleType values_;
      mutable JacobianRangeTupleType jacobians_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH
