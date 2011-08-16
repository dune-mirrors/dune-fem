#ifndef DUNE_LOCALFUNCTIONADAPTER_HH
#define DUNE_LOCALFUNCTIONADAPTER_HH

#include <set>

//- local includes 
#include <dune/fem/function/common/discretefunction.hh>

namespace Dune
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
  template <class EvalImp>
  class LocalFunctionAdapter;


  //! traits of DiscreteFunctionAdapter 
  template <class EvalImp>
  struct LocalFunctionAdapterTraits 
  {
    typedef typename EvalImp :: FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;
    typedef typename EvalImp :: GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    //! type of iterator 
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType; 
    //! type of IndexSet 
    typedef typename GridPartType :: IndexSetType IndexSetType; 
    typedef DiscreteFunctionSpaceAdapter<FunctionSpaceType,GridPartType>
            DiscreteFunctionSpaceType;

    typedef LocalFunctionAdapter<EvalImp> DiscreteFunctionType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  };

  //! identifier to local function has initialize feature 
  struct LocalFunctionAdapterHasInitialize {} ;

  /** \brief LocalFunctionAdapter wrapped a class with a local evaluate method
   *         into a grid function. 
   *
   *  The class takes one template argument 
   *  EvalImp which holds the evaluate method for the local function:
   *    template< class PointType >
   *    EvalImp::evaluate(const PointType& x,RangeType& val)
   *  Required type in EvalImp are:
   *  FunctionSpaceType which is derived from the functionspace interface provides
   *                    and provides the RangeType 
   *  GridPartType providing the EntityType
   *  An instance of the EvalImp class is passed to the constructor of the
   *  wrapper and the entity to process is passed to a method init on EvalImp.
   */
  template< class LocalFunctionImpl >
  class LocalFunctionAdapter
    : public Fem::Function< typename LocalFunctionImpl::FunctionSpaceType, LocalFunctionAdapter< LocalFunctionImpl > >
    , public HasLocalFunction
  {
    typedef LocalFunctionAdapter< LocalFunctionImpl > ThisType;
    typedef Fem::Function< typename LocalFunctionImpl::FunctionSpaceType, ThisType > BaseType;

  public:  
    typedef ThisType  DiscreteFunctionType ;

    //! Evaluate class
    typedef LocalFunctionImpl  LocalFunctionImplType ;

    //! type of function 
    typedef typename BaseType::FunctionType FunctionType;

    //! traits class
    typedef LocalFunctionAdapterTraits< LocalFunctionImplType > Traits;

    //! type of grid part 
    typedef typename LocalFunctionImplType::GridPartType GridPartType;
    
    //! type of discrete function
    typedef DiscreteFunctionSpaceAdapter< typename LocalFunctionImplType::FunctionSpaceType,
                                          GridPartType >
      DiscreteFunctionSpaceType;

    //! type of grid 
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    //! domain type (from function space)
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType ;
    //! range type (from function space)
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType ;
    //! domain type (from function space)
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType ;
    //! range type (from function space)
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType ;
    //! jacobian type (from function space)
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType 
            JacobianRangeType;
    //! type of codim 0 entity
    typedef typename GridType :: template Codim<0> :: Entity EntityType; 

  private:
    struct LocalFunction;
  public:
    //! type of local function to export 
    typedef LocalFunction LocalFunctionType; 

  protected:
    //! set of created local functions  
    typedef std::set< LocalFunctionType* > LocalFunctionListType; 

    enum { localFunctionHasInitialize = 
      Conversion< LocalFunctionImplType, LocalFunctionAdapterHasInitialize > :: exists };

    template <class ArgumentType, bool hasInit > 
    struct LocalFunctionInitializer 
    {
      static void init( const ArgumentType& , LocalFunctionListType& ) 
      {
      }
    };

    template <class ArgumentType> 
    struct LocalFunctionInitializer< ArgumentType, true >
    {
      static void init( const ArgumentType& arg, LocalFunctionListType& lfList) 
      {
        typedef typename LocalFunctionListType :: iterator iterator ;
        const iterator endit = lfList.end();
        for( iterator it = lfList.begin(); it != endit; ++it ) 
        {
          arg.initialize( (*it) );
        }
      }
    };

    // interface class for local function init 
    struct ArgumentIF 
    {
      virtual void initialize( LocalFunctionType* lf ) const = 0;
    };

    // storage of argument reference to init local functions 
    template <class ArgType> 
    struct ArgumentInitializer : public ArgumentIF 
    {
      // store arg here, this is a tuple of discrete functions 
      // that has to be copied 
      const ArgType arg_;
      // constructor storing argument  
      ArgumentInitializer( const ArgType& arg ) : arg_( arg ) {}
      virtual void initialize( LocalFunctionType* lf ) const 
      {
        lf->initialize( arg_ );
      }
    };

  public:  
    //! constructer taking instance of EvalImp class 
    LocalFunctionAdapter ( const std::string &name,
                           LocalFunctionImplType &eval,
                           const GridPartType &gridPart,
                           unsigned int order = DiscreteFunctionSpaceType::polynomialOrder )
    : space_( gridPart, order ),
      localFunctionImpl_( eval ),
      lfList_(),
      argInitializer_( 0 ),
      name_( name )
    {}

    // reference to function this local belongs to
    LocalFunctionAdapter( const ThisType &other )
    : space_( other.space_ ),
      localFunctionImpl_( other.localFunctionImpl_ ),
      lfList_(),
      argInitializer_( 0 ),
      name_( other.name_ )
    {}

    ~LocalFunctionAdapter() 
    {
      delete argInitializer_ ;
    }

    //! evaluate function on local coordinate local 
    void evaluate(const DomainType& global, RangeType& result) const 
    {
      DUNE_THROW( NotImplemented, "LocalFunctionAdapter::evaluate is not implemented." );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */ 
    const LocalFunctionType localFunction( const EntityType &entity ) const 
    {
      return LocalFunctionType( entity, *this );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) */ 
    LocalFunctionType localFunction( const EntityType &entity )
    {
      return LocalFunctionType( entity, *this );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::name */
    const std::string &name() const
    {
      return name_;
    }

    const DiscreteFunctionSpaceType &space () const
    {
      return space_;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterfaceType &g) */
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
    inline DiscreteFunctionType &operator*= ( const RangeFieldType &scalar )
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
    inline DiscreteFunctionType &operator/= ( const RangeFieldType &scalar )
    { 
      DUNE_THROW( NotImplemented, "LocalFunctionAdapter::operator /= is not implemented." );
      return *this; 
    }

    //! initialize local function with argument (see insertfunctionpass.hh)
    template <class ArgumentType>
    void initialize( const ArgumentType& arg )
    {
      if( localFunctionHasInitialize )
      {
        delete argInitializer_ ;
        // makes a copy of arg, which is a tuple of discrete functions 
        argInitializer_ = new ArgumentInitializer< ArgumentType > ( arg );
        LocalFunctionInitializer< ArgumentIF, localFunctionHasInitialize > :: init( *argInitializer_, lfList_ );
      }
      else
      {
        DUNE_THROW(NotImplemented,"LocalFunctionAdapter::initialize is not implemented");
      }
    }

    //! add LocalFunction to list of local functions 
    void registerLocalFunction( LocalFunctionType* lf ) const 
    {
      if( localFunctionHasInitialize )
      {
        if( argInitializer_ ) 
          argInitializer_->initialize( lf );
        lfList_.insert( lf );
      }
    }

    //! remove LocalFunction to list of local functions 
    void deleteLocalFunction( LocalFunctionType* lf ) const 
    {
      if( localFunctionHasInitialize )
      {
        lfList_.erase( lf );
      }
    }

  protected:    
    DiscreteFunctionSpaceType space_; 
    LocalFunctionImplType& localFunctionImpl_;
    mutable LocalFunctionListType lfList_;
    const ArgumentIF* argInitializer_ ;
    const std::string name_;
  };


  template< class LocalFunctionImpl >
  struct LocalFunctionAdapter< LocalFunctionImpl >::LocalFunction
  {
    //! domain type (from function space)
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType ;
    //! range type (from function space)
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType ;
    //! domain type (from function space)
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType ;
    //! range type (from function space)
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType ;
    //! jacobian type (from function space)
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! hessian type (from function space)
    typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;

    //! type of local function implementation 
    typedef LocalFunctionImpl LocalFunctionImplType;

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
    LocalFunction(const EntityType& en, const ThisType& a)
    : adapter_( a ),
      localFunctionImpl_(a.localFunctionImpl_)
    {
      // add local function to list
      adapter_.registerLocalFunction( this );
      localFunctionImpl_.init(en);
    }

    //! constructor 
    LocalFunction(const ThisType& a)
    : adapter_( a ),
      localFunctionImpl_(a.localFunctionImpl_) 
    {
      // add local function to list
      adapter_.registerLocalFunction( this );
    }

    //! copy constructor 
    LocalFunction(const LocalFunction& org) 
    : adapter_( org.adapter_ ),
      localFunctionImpl_(org.localFunctionImpl_) 
    {
      // add local function to list
      adapter_.registerLocalFunction( this );
    }

    //! destructor 
    ~LocalFunction() 
    {
      // remove local function from list
      adapter_.deleteLocalFunction( this );
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
    void evaluateQuadrature( const QuadratureType &quad, 
                             VectorType &result ) const
    {
      evaluateQuadrature( quad, result, result[ 0 ] );
    }

    //! init local function
    void init(const EntityType& en)
    {
      localFunctionImpl_.init(en);
    } 

    template <class ArgumentType>
    void initialize( const ArgumentType& arg )
    {
      localFunctionImpl_.initialize( arg );
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
    const ThisType& adapter_;
    typedef typename LocalFuncType<0, 
       localFunctionHasInitialize> :: Type  LocalFuncStorageType;
    LocalFuncStorageType  localFunctionImpl_;
  };

} // end namespace Dune

//@}

#endif // #ifndef DUNE_DISCRETEFUNCTIONADAPTER_HH
