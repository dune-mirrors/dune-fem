#ifndef DUNE_FEM_LOCALFUNCTIONADAPTER_HH
#define DUNE_FEM_LOCALFUNCTIONADAPTER_HH

#include <functional>
#include <memory>
#include <set>
#include <type_traits>

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/memory.hh>

#include <dune/fem/function/common/discretefunction.hh>

namespace Dune
{

  namespace Fem
  {

    /** @addtogroup DiscreteFunctionAdapter

        Similar to DiscreteFunctionAdapter but here we provide a LocalFunction
        implementation to plug this into an operator taking
        \ref DiscreteFunctionInterface "discrete functions",
        i.e., expecting \ref LocalFunction "local functions"
        a wrapper can be applied to the analytical function
        instance.
        @{
    */
    template< class LocalFunctionImpl >
    class LocalFunctionAdapter;

    template< class LocalFunctionImpl >
    class LocalFunctionAdapterLocalFunction;

    //! identifier to local function has initialize feature
    struct LocalFunctionAdapterHasInitialize {} ;



    //! traits of DiscreteFunctionAdapter
    template< class LocalFunctionImpl >
    struct LocalFunctionAdapterTraits
    {
      typedef typename LocalFunctionImpl::FunctionSpaceType FunctionSpaceType;
      typedef typename LocalFunctionImpl::GridPartType GridPartType;

      static const bool localFunctionHasInitialize = std::is_convertible< LocalFunctionImpl, LocalFunctionAdapterHasInitialize >::value;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      //! type of iterator
      typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
      //! type of IndexSet
      typedef typename GridPartType :: IndexSetType IndexSetType;

      typedef DiscreteFunctionSpaceAdapter< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;

      typedef LocalFunctionAdapter< LocalFunctionImpl > DiscreteFunctionType;
      typedef LocalFunctionAdapterLocalFunction< LocalFunctionImpl > LocalFunctionType;
    };



    /** \brief LocalFunctionAdapter wrapped a class with a local evaluate method
     *         into a grid function.
     *
     *  The class takes one template argument LocalFunctionImpl which holds the
     *  method evaluate(...) to evaluate the local function
     *
     *    template<class PointType>
     *    LocalFunctionImpl::evaluate(const PointType& x,RangeType& val)
     *
     *  and a method init(...)
     *
     *    LocalFunctionImpl::init(const EntityType& entity)
     *
     *  to set the entity.
     *
     *  It is important to know that the point x it is not necessary of type DomainType.
     *  More precisely, if the evaluate(...) is used with a caching quadrature point the
     *  type is different. Indeed floating point coordinates are not very well suited to
     *  address the cache therefore quadrature[i] return a QuadraturePointWrapper
     *  (which simply stores a reference to the quadrature and the index i).
     *
     *  In order to be sure that the point x is of type DomainType, you can use the function
     *  coordinate(x) which can be also called with a DomainType.
     *
     *  Therefore, the local implementation should be something like
     *
     *    template<class PointType>
     *    LocalFunctionImpl::evaluate(const PointType& x,RangeType& val)
     *    {
     *      const DomainType xDomain(coordiante(x));
     *      // do stuff with xDomain
     *    }
     *
     *  to avoid type conflicts.
     *
     *  Required type in LocalFunctionImpl are:
     *
     *  FunctionSpaceType
     *  GridPartType
     *  EntityType
     *  DomainType
     *  RangeType
     *
     *  An instance of the LocalFunctionImpl class is passed to the constructor.
     *
     *  In order to adapt a lambda or a plain C++ function, you can directly use
     *  the LocalAnalyticalFunctionBinder which provides all the necessary
     *  types and methods.
     */
    template< class LocalFunctionImpl >
    class LocalFunctionAdapter
    : public Function< typename LocalFunctionImpl::FunctionSpaceType, LocalFunctionAdapter< LocalFunctionImpl > >,
      public HasLocalFunction
    {
      typedef LocalFunctionAdapter< LocalFunctionImpl > ThisType;
      typedef Function< typename LocalFunctionImpl::FunctionSpaceType, ThisType > BaseType;

      friend class LocalFunctionAdapterLocalFunction< LocalFunctionImpl >;

    public:
      typedef ThisType  DiscreteFunctionType;

      //! Evaluate class
      typedef LocalFunctionImpl LocalFunctionImplType;

      //! type of function
      typedef typename BaseType::FunctionType FunctionType;

      //! traits class
      typedef LocalFunctionAdapterTraits< LocalFunctionImplType > Traits;

      //! type of grid part
      typedef typename Traits::GridPartType GridPartType;

      //! type of discrete function space
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      //! type of grid
      typedef typename Traits::GridType GridType;
      //! domain type
      typedef typename Traits::DomainFieldType DomainFieldType ;
      //! range type
      typedef typename Traits::RangeFieldType RangeFieldType ;
      //! domain type
      typedef typename Traits::DomainType DomainType;
      //! range type
      typedef typename Traits::RangeType RangeType;
      //! jacobian type
      typedef typename Traits::JacobianRangeType JacobianRangeType;
      //! type of codim 0 entity
      typedef typename Traits::EntityType EntityType;

      //! type of local function to export
      typedef typename Traits::LocalFunctionType LocalFunctionType;

    protected:
      // interface class for local function init
      struct ArgumentIF
      {
        virtual void initialize( LocalFunctionType* ) const = 0;

        virtual ~ArgumentIF ()
        {}
      };

      // storage of argument reference to init local functions
      template <class ArgType>
      struct ArgumentInitializer : public ArgumentIF
      {
        // store arg here, this is a tuple of discrete functions
        // that has to be copied
        const ArgType arg_;
        const double time_;

        // constructor storing argument
        ArgumentInitializer( const ArgType& arg, double time )
        : arg_( arg ), time_( time )
        {}

        ~ArgumentInitializer ()
        {}

        virtual void initialize( LocalFunctionType* lf ) const
        {
          lf->initialize( arg_, time_ );
        }
      };

    public:
      //! constructer taking instance of EvalImp class
      LocalFunctionAdapter ( const std::string &name,
                             LocalFunctionImplType &localFunctionImpl,
                             const GridPartType &gridPart,
                             unsigned int order = DiscreteFunctionSpaceType::polynomialOrder )
      : space_( gridPart, order ),
        localFunctionImpl_( localFunctionImpl ),
        lfList_(),
        argInitializer_(),
        name_( name ),
        order_( order )
      {}

      // reference to function this local belongs to
      LocalFunctionAdapter( const ThisType &other )
      : space_( other.space_ ),
        localFunctionImpl_( other.localFunctionImpl_ ),
        lfList_(),
        argInitializer_(),
        name_( other.name_ ),
        order_( other.order_ )
      {}

      //! return the order of the space
      unsigned int order() const
      {
        return order_;
      }

      //! return local function implementation
      const LocalFunctionImplType& localFunctionImpl() const
      {
        return localFunctionImpl_;
      }

      //! return local function implementation
      LocalFunctionImplType& localFunctionImpl()
      {
        return localFunctionImpl_;
      }

      //! evaluate function on local coordinate local
      void evaluate(const DomainType& global, RangeType& result) const
      {
        DUNE_THROW( NotImplemented, "LocalFunctionAdapter::evaluate is not implemented." );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      LocalFunctionType localFunction( const EntityType &entity )
      {
        return LocalFunctionType( entity, *this );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      const LocalFunctionType localFunction( const EntityType &entity ) const
      {
        return LocalFunctionType( entity, *this );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::name */
      const std::string &name() const
      {
        return name_;
      }

      const DiscreteFunctionSpaceType &space () const
      {
        return space_;
      }

      const GridPartType &gridPart () const
      {
        return space().gridPart();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterfaceType &g) */
      template <class DFType>
      DiscreteFunctionType &operator+= ( const DFType &g )
      {
        DUNE_THROW( NotImplemented, "LocalFunctionAdapter::operator += is not implemented." );
        return *this;
      }

      /** \brief substract all degrees of freedom from given discrete function using the dof iterators
          \param[in] g discrete function which is substracted from this discrete function
          \return reference to this (i.e. *this)
      */
      template <class DFType>
      DiscreteFunctionType& operator -= (const DFType& g)
      {
        DUNE_THROW( NotImplemented, "LocalFunctionAdapter::operator -= is not implemented." );
        return *this;
      }

      /** \brief multiply all DoFs with a scalar factor
       *
       *  \param[in]  scalar  factor to multiply DoFs with
       *
       *  \returns reference to this discrete function (i.e. *this)
       */
      DiscreteFunctionType &operator*= ( const RangeFieldType &scalar )
      {
        DUNE_THROW( NotImplemented, "LocalFunctionAdapter::operator *= is not implemented." );
        return *this;
      }

      /** \brief devide all DoFs by a scalar factor
       *
       *  \param[in]  scalar  factor with which all dofs are devided
       *
       *  \returns reference to this discrete function (i.e. *this)
       */
      DiscreteFunctionType &operator/= ( const RangeFieldType &scalar )
      {
        DUNE_THROW( NotImplemented, "LocalFunctionAdapter::operator /= is not implemented." );
        return *this;
      }

      //! initialize local function with argument
      template< class ArgumentType >
      void initialize( const ArgumentType &arg, double time )
      {
        constexpr bool hasInit = Traits::localFunctionHasInitialize;
        if( hasInit)
        {
          // makes a copy of arg, which is a tuple of discrete functions
          argInitializer_ = Std::make_unique< ArgumentInitializer< ArgumentType > >( arg, time );
          for( auto& localFunctionPtr : lfList_ )
            Hybrid::ifElse( hasInit, [ & ]( auto&& ){ argInitializer_->initialize( localFunctionPtr ); } );
        }
        else
          DUNE_THROW(NotImplemented,"LocalFunctionAdapter::initialize is not implemented");
      }

      //! add LocalFunction to list of local functions
      void registerLocalFunction( LocalFunctionType* lf ) const
      {
        if( Traits::localFunctionHasInitialize )
        {
          if( argInitializer_ != nullptr )
            argInitializer_->initialize( lf );
          lfList_.insert( lf );
        }
      }

      //! remove LocalFunction to list of local functions
      void deleteLocalFunction( LocalFunctionType* lf ) const
      {
        if( Traits::localFunctionHasInitialize )
          lfList_.erase( lf );
      }

    protected:
      DiscreteFunctionSpaceType space_;
      LocalFunctionImplType &localFunctionImpl_;
      mutable std::set< LocalFunctionType * > lfList_;
      std::unique_ptr< ArgumentIF > argInitializer_ ;
      const std::string name_;
      const unsigned int order_;
    };



    template< class LocalFunctionImpl >
    class LocalFunctionAdapterLocalFunction
    {
      typedef LocalFunctionAdapterLocalFunction< LocalFunctionImpl > ThisType;

    public:
      //! type of local function implementation
      typedef LocalFunctionImpl LocalFunctionImplType;

      //! type of the traits class
      typedef LocalFunctionAdapterTraits< LocalFunctionImplType > Traits;

      //! domain type
      typedef typename Traits::DomainFieldType DomainFieldType;
      //! range type
      typedef typename Traits::RangeFieldType RangeFieldType;
      //! domain type
      typedef typename Traits::DomainType DomainType;
      //! range type
      typedef typename Traits::RangeType RangeType;
      //! jacobian type
      typedef typename Traits::JacobianRangeType JacobianRangeType;
      //! hessian type
      typedef typename Traits::HessianRangeType HessianRangeType;

      typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
      typedef typename Traits::EntityType EntityType;

      // default is reference
      template <int, bool hasInit >
      struct LocalFuncType
      {
        typedef LocalFunctionImplType&  Type;
      };

      // non default is object
      template <int dummy >
      struct LocalFuncType<dummy, true>
      {
        typedef LocalFunctionImplType Type;
      };

      //! constructor initializing local function
      LocalFunctionAdapterLocalFunction ( const EntityType &entity, const DiscreteFunctionType &adapter )
      : adapter_( adapter ), localFunctionImpl_( adapter.localFunctionImpl_ )
      {
        // add local function to list
        adapter_.registerLocalFunction( this );
        localFunctionImpl_.init( entity );
      }

      //! constructor
      explicit LocalFunctionAdapterLocalFunction ( const DiscreteFunctionType &adapter )
      : adapter_( adapter ), localFunctionImpl_( adapter.localFunctionImpl_ )
      {
        // add local function to list
        adapter_.registerLocalFunction( this );
      }

      //! copy constructor
      LocalFunctionAdapterLocalFunction ( const ThisType &other )
      : adapter_( other.adapter_ ), localFunctionImpl_( other.localFunctionImpl_ )
      {
        // add local function to list
        adapter_.registerLocalFunction( this );
      }

      //! destructor
      ~LocalFunctionAdapterLocalFunction ()
      {
        // remove local function from list
        adapter_.deleteLocalFunction( this );
      }

      //! return order of the space
      unsigned int order() const
      {
        return adapter_.order();
      }

      //! evaluate local function
      template< class PointType >
      void evaluate ( const PointType &x, RangeType &ret ) const
      {
        localFunctionImpl_.evaluate(x,ret);
      }

      //! jacobian of local function
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        localFunctionImpl_.jacobian( x, ret );
      }

      // hessian of local function
      template< class PointType >
      void hessian ( const PointType &x, HessianRangeType &ret ) const
      {
        localFunctionImpl_.hessian( x, ret );
      }

      template< class QuadratureType, class VectorType  >
      void evaluateQuadrature( const QuadratureType &quad, VectorType &result ) const
      {
        evaluateQuadrature( quad, result, result[ 0 ] );
      }

      //! init local function
      void init(const EntityType& en)
      {
        localFunctionImpl_.init(en);
        entity_=&en;
      }

      template <class ArgumentType>
      void initialize ( const ArgumentType& arg, double time )
      {
        localFunctionImpl_.initialize( arg, time );
      }

      //! get entity
      const EntityType& entity() const
      {
        assert( entity_ );
        return *entity_;
      }

    protected:
      template< class QuadratureType, class VectorType  >
      void evaluateQuadrature( const QuadratureType &quad,
          VectorType &result, const RangeType& ) const
      {
        const size_t quadNop = quad.nop();
        for(size_t i = 0; i<quadNop; ++i)
          evaluate( quad[ i ], result[ i ] );
      }

      template< class QuadratureType, class VectorType  >
      void evaluateQuadrature( const QuadratureType &quad,
          VectorType &result, const JacobianRangeType& ) const
      {
        const size_t quadNop = quad.nop();
        for(size_t i = 0; i<quadNop; ++i)
          jacobian( quad[ i ], result[ i ] );
      }

    protected:
      const DiscreteFunctionType &adapter_;
      typedef typename LocalFuncType< 0, Traits::localFunctionHasInitialize >::Type LocalFuncStorageType;
      LocalFuncStorageType localFunctionImpl_;
      EntityType const* entity_;
    };



    /** \brief LocalAnalyticalFunctionBinder binds a C++ local analytical function (and also its Jacobian
     *  and Hessian) to an object which provides all the methods and types needed by the LocalFunctionAdapter.
     *  It stores a copy as a std::function.
     *
     *  Therefore, in order to transform the function
     *
     *    RangeType f(const DomainType& x,double t,const EntityType& entity)
     *    {
     *      // do stuff
     *    }
     *
     *  into a grid function, it is sufficient to pass it to the LocalAnalyticalFucntionBinder
     *
     *    typedef LocalAnalyticalFunctionBinder<DiscreteFunctionSpaceType> LocalAnalyticalFunctionType;
     *    LocalAnalyticalFunctionType localAnalyticalFunction(f);
     *
     *  and create the LocalFunctionAdapter
     *
     *    typedef LocalFunctionAdapter<LocalAnalyticalFunctionType> AdaptedFunctionType;
     *    AdaptedFunctionType fAdapted("adapted function",localAnalyticalFunction,gridPart);
     */
    template<class DiscreteFunctionSpaceImpl>
    class LocalAnalyticalFunctionBinder:public LocalFunctionAdapterHasInitialize
    {
    public:
      typedef DiscreteFunctionSpaceImpl DiscreteFunctionSpaceType;
      typedef LocalAnalyticalFunctionBinder<DiscreteFunctionSpaceType> ThisType;

      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef std::function<RangeType(const DomainType&,double,const EntityType&)> AnalyticalFunctionType;
      typedef std::function<JacobianRangeType(const DomainType&,double,const EntityType&)> AnalyticalJacobianType;
      typedef std::function<HessianRangeType(const DomainType&,double,const EntityType&)> AnalyticalHessianType;

      //! constructor
      LocalAnalyticalFunctionBinder(const AnalyticalFunctionType& f=[](const auto& ,auto ,const auto& ){return RangeType(0.0);},
                                    const AnalyticalJacobianType& j=[](const auto& ,auto ,const auto& ){return JacobianRangeType(0.0);},
                                    const AnalyticalHessianType& h=[](const auto& ,auto ,const auto& ){return HessianRangeType(0.0);},
                                    double t=0.0):
        f_(f),j_(j),h_(h),t_(t)
      {}

      LocalAnalyticalFunctionBinder(const ThisType& )=default;
      LocalAnalyticalFunctionBinder(ThisType&& )=default;
      ThisType& operator=(const ThisType& )=default;
      ThisType& operator=(ThisType&& )=default;

      //! set local function
      void setFunction(const AnalyticalFunctionType& f)
      {
        f_=f;
      }

      //! set jacobian local function
      void setJacobian(const AnalyticalJacobianType& j)
      {
        j_=j;
      }

      //! set hessian local function
      void setHessian(const AnalyticalHessianType& h)
      {
        h_=h;
      }

      //! evaluate local function
      template<class PointType>
      void evaluate(const PointType& x,RangeType& ret) const
      {
        ret=f_(entity().geometry().global(coordinate(x)),t_,entity());
      }

      //! evaluate jacobian local function
      template<class PointType>
      void jacobian(const PointType& x,JacobianRangeType& ret) const
      {
        ret=j_(entity().geometry().global(coordinate(x)),t_,entity());
      }

      //! evaluate hessian local function
      template<class PointType>
      void hessian(const PointType& x,HessianRangeType& ret ) const
      {
        ret=h_(entity().geometry().global(coordinate(x)),t_,entity());
      }

      //! initialize entity
      void init(const EntityType& entity)
      {
        entity_=&entity;
      }

      //! initialize time
      template<typename Arg>
      void initialize(Arg&& ,double time)
      {
        t_=time;
      }

      //! get entity
      const EntityType& entity() const
      {
        assert(entity_);
        return *entity_;
      }

      //! get time
      double time() const
      {
        return t_;
      }

    private:
      EntityType const* entity_;
      AnalyticalFunctionType f_;
      AnalyticalJacobianType j_;
      AnalyticalHessianType h_;
      double t_;
    };

  } // namespace Fem

} // namespace Dune

//@}

#endif // #ifndef DUNE_FEM_LOCALFUNCTIONADAPTER_HH
