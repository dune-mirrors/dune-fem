#ifndef DUNE_DISCRETEFUNCTIONADAPTER_HH
#define DUNE_DISCRETEFUNCTIONADAPTER_HH

//- system includes 
#include <string>

//- local includes 
#include <dune/fem/function/common/discretefunction.hh>

namespace Dune{

  /** 
      @addtogroup DiscreteFunctionAdapter

      To plug an \ref Function "analytical function" 
      into a operator taking 
      \ref DiscreteFunctionInterface "discrete functions",
      i.e., expecting \ref LocalFunction "local functions"
      a wrapper can be applied to the analytical function
      instance.
      The resulting class is still a \ref Function "Function"
      but with the property \ref HasLocalFunction "\" has local function \"" added.
      
      @{
  */

  //- forward declaration 
  template <class FunctionImp, class GridPartImp>
  class DiscreteFunctionAdapter;

  //! traits of DiscreteFunctionAdapter 
  template <class FunctionImp, class GridPartImp> 
  struct DiscreteFunctionAdapterTraits 
  {
    // make sure we obtain only the function space
    typedef typename FunctionImp :: FunctionSpaceType :: FunctionSpaceType
      FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef DiscreteFunctionSpaceAdapter<FunctionSpaceType,GridPartImp> DiscreteFunctionSpaceType;

    typedef GridPartImp GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    //! type of iterator 
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType; 
    //! type of IndexSet 
    typedef typename GridPartType :: IndexSetType IndexSetType; 

    typedef DiscreteFunctionAdapter<FunctionImp,GridPartImp> DiscreteFunctionType;
  };

  /** \brief DiscreteFunctionAdapter provides local functions for a Function. 
   */
  template <class FunctionImp, class GridPartImp>
  class DiscreteFunctionAdapter 
  : public HasLocalFunction , 
    public Function<DiscreteFunctionSpaceAdapter<typename FunctionImp :: FunctionSpaceType,GridPartImp>, 
                    DiscreteFunctionAdapter<FunctionImp,GridPartImp> >
  {
  public:  
    //! type of function 
    typedef FunctionImp FunctionType;

    //! type of grid part 
    typedef GridPartImp GridPartType;
                       
    //! type of traits 
    typedef DiscreteFunctionAdapterTraits< FunctionType, GridPartType > Traits;

    //! type of discrete function space 
    typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  private:
    typedef DiscreteFunctionAdapter< FunctionType, GridPartType > ThisType;
    typedef Function< DiscreteFunctionSpaceType, ThisType > BaseType;

    // Make sure the function is not a discrete functon
    CompileTimeChecker < !(Conversion< FunctionType, HasLocalFunction > :: exists) >
      __FunctionType_May_Not_Be_a_DiscreteFunctionType__;

  public:
    // type of discrete function space
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType; 

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
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    //! type of codim 0 entity
    typedef typename GridType :: template Codim<0> :: Entity EntityType; 

    private:
    class LocalFunctionStorage;
    class LocalFunction
    {
      //! type of geometry 
      typedef typename EntityType :: Geometry GeometryImp;
    public:  
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

      //! constructor initializing local function 
      LocalFunction(const EntityType& en, const ThisType& a)
        : function_(a.function_) 
        , geometry_(&(en.geometry())) 
      {}

      LocalFunction(const ThisType& a)
        : function_(a.function_) 
        , geometry_(0) 
      {}

      //! copy constructor 
      LocalFunction(const LocalFunction& org) 
        : function_(org.function_) 
        , geometry_(org.geometry_)  
      {}
      
      LocalFunction(LocalFunctionStorage& storage) 
        : function_(storage.function().function_)
        , geometry_(0) 
      {}

      //! evaluate local function 
      template< class PointType >
      void evaluate ( const PointType &x, RangeType &ret ) const
      {
        DomainType global = geometry_->global( coordinate( x ) );
        function_.evaluate( global, ret );
      }

      //! jacobian of local function 
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        DomainType global = geometry_->global( coordinate( x ) );
        function_.jacobian( global, ret );
      }

      //! init local function
      void init(const EntityType& en) 
      {
        geometry_ = &(en.geometry());
      } 

    private:
      const FunctionType& function_;
      const GeometryImp* geometry_;
    };


    public:
    //! type of local function to export 
    typedef LocalFunction LocalFunctionType; 
    typedef LocalFunctionStorage LocalFunctionStorageType;

    // reference to function this local belongs to
    inline DiscreteFunctionAdapter
      ( const std :: string &name,
        const FunctionType &f,
        const GridPartType &gridPart,
        unsigned int order = DiscreteFunctionSpaceType :: polynomialOrder )
    : BaseType(space_),
      space_( gridPart, order ),
      localFunctionStorage_( *this ),
      function_( f ),
      name_( name )
    {
    }

    // reference to function this local belongs to
    DiscreteFunctionAdapter( const ThisType &other ) 
    : BaseType( other ),
      space_( other.space_ ),
      localFunctionStorage_( *this ),
      function_( other.function_ ),
      name_( other.name_ )
    {
    }

    //! evaluate function on local coordinate local 
    void evaluate(const DomainType& global, RangeType& result) const 
    {
      function_.evaluate(global,result);  
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */ 
    const LocalFunctionType localFunction( const EntityType &entity ) const 
    {
      return localFunctionStorage().localFunction( entity );
      //return LocalFunctionType( entity, *this );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) */ 
    LocalFunctionType localFunction( const EntityType &entity )
    {
      return localFunctionStorage().localFunction( entity );
      //return LocalFunctionType( entity, *this );
    }

    inline LocalFunctionStorageType &localFunctionStorage () const
    {
      return localFunctionStorage_;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::name */
    inline const std :: string &name() const
    {
      return name_;
    }

  protected:    
    DiscreteFunctionSpaceType space_;
    mutable LocalFunctionStorageType localFunctionStorage_;
    //! reference to function 
    const FunctionType& function_; 
    
    const std::string name_;
  };



  template< class Function, class GridPart >
  class DiscreteFunctionAdapter< Function, GridPart > :: LocalFunctionStorage
  {
    typedef LocalFunctionStorage ThisType;
    typedef DiscreteFunctionAdapter< Function, GridPart > DiscreteFunctionType;

  public:
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

  private:
    DiscreteFunctionType &discreteFunction_;

  public:
    inline explicit
    LocalFunctionStorage ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction )
    {}

  private:
    LocalFunctionStorage ( const ThisType & );
    ThisType operator= ( const ThisType & );

  public:
    inline LocalFunctionType localFunction ()
    {
      return LocalFunctionType( discreteFunction_ );
    }

    template< class Entity >
    inline const LocalFunctionType localFunction ( const Entity &entity ) const
    {
      return LocalFunctionType( entity, discreteFunction_ );
    }

    template< class Entity >
    inline LocalFunctionType localFunction ( const Entity &entity )
    {
      return LocalFunctionType( entity, discreteFunction_ );
    }

    DiscreteFunctionType& function() {
      return discreteFunction_;
    }
  };

  

  namespace {
    template <class FunctionImp,class GridPartType,bool>
    struct ConvertDFTypeHelper;
    template <class FunctionImp,class GridPartType>
    struct ConvertDFTypeHelper<FunctionImp,GridPartType,true> 
    {
      typedef ConvertDFTypeHelper<FunctionImp,GridPartType,true> 
        ThisType;
      enum {compatible =
        Conversion<GridPartType,typename
          FunctionImp::DiscreteFunctionSpaceType::GridPartType>::exists};
      typedef FunctionImp FunctionType;
      typedef typename FunctionType::DiscreteFunctionSpaceType DFSType;
      ConvertDFTypeHelper(const std::string& name,const FunctionImp& func,const GridPartType& gp) :
        func_(func) {}
      ConvertDFTypeHelper(const ConvertDFTypeHelper& other) :
        func_(other.func_) {}
      const FunctionType& function() const {
        return func_;
      }
      const DFSType& space() const {
        return func_.space();
      }
      private:
      const FunctionImp& func_;
    };
    template <class FunctionImp,class GridPartType>
    struct ConvertDFTypeHelper<FunctionImp,GridPartType,false> :
      DiscreteFunctionAdapter<FunctionImp,GridPartType>
    {
      typedef ConvertDFTypeHelper<FunctionImp,GridPartType,false>
        ThisType;
      typedef DiscreteFunctionAdapter<FunctionImp,GridPartType> BaseType;
      typedef BaseType FunctionType;
      typedef typename FunctionType::DiscreteFunctionSpaceType DFSType;
      ConvertDFTypeHelper(const std::string& name,const FunctionImp& func,const GridPartType& gp) :
        BaseType(name,func,gp) {}
      ConvertDFTypeHelper(const ConvertDFTypeHelper& other) :
        BaseType(other) {}
      const FunctionType& function() const {
        return *this;
      }
      const DFSType& space() const {
        return BaseType::space();
      }
    };
  }
  template <class FunctionImp,class GridPartImp>
  class ConvertToGridFunction : 
    public Function<typename
           ConvertDFTypeHelper<FunctionImp,GridPartImp,
           Conversion< FunctionImp, HasLocalFunction > :: exists
           >::DFSType,ConvertToGridFunction<FunctionImp,GridPartImp>
           >,
    public HasLocalFunction
  {
  public:
    typedef FunctionImp FunctionType;
    typedef GridPartImp GridPartType;
  private:
    enum { hasLocalFunction = Conversion< FunctionType, HasLocalFunction > :: exists  };
    typedef ConvertToGridFunction<FunctionType,GridPartImp>
      ThisType;
    typedef ConvertDFTypeHelper<FunctionImp,GridPartImp,hasLocalFunction>
      Helper;
    typedef Function<typename Helper::DFSType,ThisType> BaseType;
    typedef typename Helper::FunctionType ConvertedType;
  public:
    //! type of discrete function space 
    typedef typename ConvertedType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    // type of discrete function space
    typedef typename ConvertedType :: FunctionSpaceType FunctionSpaceType; 

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
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    //! type of codim 0 entity
    typedef typename GridType :: template Codim<0> :: Entity EntityType; 

    //! type of local function to export 
    typedef typename ConvertedType::LocalFunctionType LocalFunctionType; 

    //! constructor
    ConvertToGridFunction(const std::string& name,
                          const FunctionImp& func,
                          const GridPartType& gridPart) :
      BaseType(helper_.space()),
      name_(name),
      helper_(name,func,gridPart) {}
    ConvertToGridFunction( const ThisType &other ) 
    : BaseType(other),
      name_(other.name_),
      helper_(other.helper_)
    {}

    //! evaluate function on local coordinate local 
    void evaluate(const DomainType& global, RangeType& result) const 
    {
      helper_.function().evaluate(global,result);  
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */ 
    const LocalFunctionType localFunction( const EntityType &entity ) const 
    {
      return helper_.function().localFunction(entity);
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) */ 
    LocalFunctionType localFunction( const EntityType &entity )
    {
      return helper_.function().localFunction(entity);
    }

    /** \copydoc Dune::DiscreteFunctionInterface::name */
    inline const std :: string &name() const
    {
      return name_;
    }
  private:
    const std::string name_;
    Helper helper_;
  };
  template <class FunctionImp,class GridPartImp>
  ConvertToGridFunction<FunctionImp,GridPartImp>
  convertToGridFunction(const std::string& name,
                        const FunctionImp& func,
                        const GridPartImp& gridPart) {
    return ConvertToGridFunction<FunctionImp,GridPartImp>
      (name,func,gridPart);
  }

  
  //! traits of DiscreteFunctionAdapter 
  template <class EvalImp>
  class LocalFunctionAdapter;
  template <class EvalImp>
  struct LocalFunctionAdapterTraits 
  {
    typedef typename EvalImp :: FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
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
  };

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
  template <class EvalImp>
  class LocalFunctionAdapter
  : public HasLocalFunction , 
    public Function<DiscreteFunctionSpaceAdapter<typename EvalImp :: FunctionSpaceType,typename EvalImp :: GridPartType>, 
                    LocalFunctionAdapter<EvalImp> >
  {
  private:
    typedef LocalFunctionAdapter<EvalImp> ThisType;
  public:  
    //! type of function 
    typedef ThisType FunctionType;

    //! traits class
    typedef LocalFunctionAdapterTraits<EvalImp> Traits;

    //! type of grid part 
    typedef typename EvalImp::GridPartType GridPartType;
    
    //! type of discrete function
    typedef DiscreteFunctionSpaceAdapter
            <typename EvalImp :: FunctionSpaceType,
	     typename EvalImp :: GridPartType>
            DiscreteFunctionSpaceType;
  private:
    typedef Function< DiscreteFunctionSpaceType, ThisType > BaseType;
  public:
    //! Evaluate class
    typedef EvalImp EvalType;
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
    class LocalFunction
    {
    public:  
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

      //! constructor initializing local function 
      LocalFunction(const EntityType& en, const ThisType& a)
        : eval_(a.eval_) {
	eval_.init(en);
      }

      LocalFunction(const ThisType& a)
        : eval_(a.eval_) 
      {}

      //! copy constructor 
      LocalFunction(const LocalFunction& org) 
        : eval_(org.eval_) 
      {}

      //! evaluate local function 
      template< class PointType >
      void evaluate ( const PointType &x, RangeType &ret ) const {
        eval_.evaluate(x,ret);
      }

      //! jacobian of local function 
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const {
        eval_.jacobian( x, ret );
      }

      //! init local function
      void init(const EntityType& en) {
        eval_.init(en);
      } 
    private:
      EvalType& eval_;
    };
    public:
    //! type of local function to export 
    typedef LocalFunction LocalFunctionType; 

    //! constructer taking instance of EvalImp class 
    inline LocalFunctionAdapter
      ( const std :: string &name,
        EvalType &eval,
        const GridPartType &gridPart,
        unsigned int order = DiscreteFunctionSpaceType :: polynomialOrder )
    : BaseType(space_),
      space_( gridPart, order ),
      eval_( eval ),
      name_( name )
    {
    }

    // reference to function this local belongs to
    LocalFunctionAdapter( const ThisType &other ) 
    : BaseType( other ),
      space_( other.space_ ),
      eval_( other.eval_ ),
      name_( other.name_ )
    {
    }

    //! evaluate function on local coordinate local 
    void evaluate(const DomainType& global, RangeType& result) const 
    {
      // DUNE_NOT_IMPLEMENTED;
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
    inline const std :: string &name() const
    {
      return name_;
    }

  private:    
    DiscreteFunctionSpaceType space_; 
    mutable EvalType& eval_;
    const std::string name_;
  };

} // end namespace Dune

//@}

#endif
