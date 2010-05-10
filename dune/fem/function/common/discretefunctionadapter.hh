#ifndef DUNE_DISCRETEFUNCTIONADAPTER_HH
#define DUNE_DISCRETEFUNCTIONADAPTER_HH

//- local includes 
#include <dune/fem/function/common/discretefunction.hh>

namespace Dune
{

  /** @addtogroup DiscreteFunctionAdapter

      To plug an \ref Fem::Function "analytical function" 
      into a operator taking 
      \ref DiscreteFunctionInterface "discrete functions",
      i.e., expecting \ref LocalFunction "local functions"
      a wrapper can be applied to the analytical function
      instance.
      The resulting class is still a \ref Fem::Function "Function"
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
    typedef typename FunctionImp::FunctionSpaceType FunctionSpaceType;

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



  // DiscreteFunctionAdapter
  // -----------------------

  /** \brief DiscreteFunctionAdapter provides local functions for a Function. 
   */
  template< class FunctionImp, class GridPartImp >
  class DiscreteFunctionAdapter 
  : public Fem::Function< typename FunctionImp::FunctionSpaceType,
                          DiscreteFunctionAdapter< FunctionImp, GridPartImp > >,
    public HasLocalFunction
  {
    typedef DiscreteFunctionAdapter< FunctionImp, GridPartImp > ThisType;
    typedef Fem::Function< typename FunctionImp::FunctionSpaceType, ThisType > BaseType;

    // Make sure the function is not a discrete functon
    dune_static_assert( !(Conversion< FunctionImp, HasLocalFunction >::exists),
                        "FunctionType may not be a discrete function type." );

  public:  
    //! type of traits 
    typedef DiscreteFunctionAdapterTraits< FunctionImp, GridPartImp > Traits;

    //! type of function 
    typedef FunctionImp FunctionType;

    //! type of grid part 
    typedef GridPartImp GridPartType;
                       
    //! type of discrete function space 
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // type of discrete function space
    typedef typename Traits::FunctionSpaceType FunctionSpaceType; 

    //! type of grid 
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

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
    typedef typename GridType::template Codim< 0 >::Entity EntityType; 

  private:
    class LocalFunction;
    class LocalFunctionStorage;

    public:
    //! type of local function to export 
    typedef LocalFunction LocalFunctionType; 
    typedef LocalFunctionStorage LocalFunctionStorageType;

    // reference to function this local belongs to
    DiscreteFunctionAdapter ( const std::string &name,
                              const FunctionType &f,
                              const GridPartType &gridPart,
                              unsigned int order = DiscreteFunctionSpaceType::polynomialOrder )
    : space_( gridPart, order ),
      localFunctionStorage_( *this ),
      function_( f ),
      name_( name )
    {}

    // reference to function this local belongs to
    DiscreteFunctionAdapter( const ThisType &other ) 
    : space_( other.space_ ),
      localFunctionStorage_( *this ),
      function_( other.function_ ),
      name_( other.name_ )
    {}

    //! evaluate function on local coordinate local 
    void evaluate(const DomainType& global, RangeType& result) const 
    {
      function_.evaluate(global,result);  
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */ 
    const LocalFunctionType localFunction( const EntityType &entity ) const 
    {
      return localFunctionStorage().localFunction( entity );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) */ 
    LocalFunctionType localFunction( const EntityType &entity )
    {
      return localFunctionStorage().localFunction( entity );
    }

    LocalFunctionStorageType &localFunctionStorage () const
    {
      return localFunctionStorage_;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::name() const */
    const std::string &name () const
    {
      return name_;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::space() const */
    const DiscreteFunctionSpaceType &space () const
    {
      return space_;
    }

  protected:    
    DiscreteFunctionSpaceType space_;
    mutable LocalFunctionStorageType localFunctionStorage_;
    const FunctionType& function_; 
    const std::string name_;
  };



  // DiscreteFunctionAdapter::LocalFunction
  // --------------------------------------

  template< class Function, class GridPart >
  class DiscreteFunctionAdapter< Function, GridPart >::LocalFunction
  {
    typedef LocalFunction ThisType;
    typedef DiscreteFunctionAdapter< Function, GridPart > DiscreteFunctionType;

    typedef typename EntityType::Geometry GeometryType;

  public:
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;
    static const int dimDomain = GridPart::GridType::dimensionworld;
    static const int dimLocal = GridPart::GridType::dimension;

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
    LocalFunction ( const EntityType &entity, const DiscreteFunctionType &df )
    : function_( &df.function_ ),
      geometry_( &entity.geometry() )
    {}

    LocalFunction ( const DiscreteFunctionType &df )
    : function_( &df.function_ ),
      geometry_( 0 )
    {}

    LocalFunction ( LocalFunctionStorage &storage )
    : function_( &storage.function().function_ ),
      geometry_( 0 )
    {}

    //! evaluate local function 
    template< class PointType >
    void evaluate ( const PointType &x, RangeType &ret ) const
    {
      DomainType global = geometry().global( coordinate( x ) );
      function().evaluate( global, ret );
    }

    //! jacobian of local function 
    template< class PointType >
    void jacobian ( const PointType &x, JacobianRangeType &ret ) const
    {
      const typename GeometryType::LocalCoordinate cx = coordinate( x );
      DomainType global = geometry().global( cx );
      function().jacobian( global, ret );

      if( dimLocal != dimDomain )
      {
        const typename GeometryType::JacobianTransposed &gjt = geometry().jacobianTransposed( cx );
        const typename GeometryType::Jacobian &gjit = geometry().jacobianInverseTransposed( cx );

        FieldVector< RangeFieldType, dimLocal > tmp;
        for( int i = 0; i < dimRange; ++i )
        {
          gjit.mtv( ret[ i ], tmp );
          gjt.mtv( tmp, ret[ i ] );
        }
      }
    }

    //! init local function
    void init ( const EntityType &entity )
    {
      geometry_ = &entity.geometry();
    } 

  private:
    const FunctionType &function () const
    {
      return *function_;
    }

    const GeometryType &geometry () const
    {
      return *geometry_;
    }

    const FunctionType *function_;
    const GeometryType *geometry_;
  };



  // DiscreteFunctionAdapter::LocalFunctionStorage
  // ---------------------------------------------

  template< class Function, class GridPart >
  class DiscreteFunctionAdapter< Function, GridPart >::LocalFunctionStorage
  {
    typedef LocalFunctionStorage ThisType;
    typedef DiscreteFunctionAdapter< Function, GridPart > DiscreteFunctionType;

  public:
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

    explicit LocalFunctionStorage ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction )
    {}

    LocalFunctionType localFunction ()
    {
      return LocalFunctionType( discreteFunction_ );
    }

    template< class Entity >
    const LocalFunctionType localFunction ( const Entity &entity ) const
    {
      return LocalFunctionType( entity, discreteFunction_ );
    }

    template< class Entity >
    LocalFunctionType localFunction ( const Entity &entity )
    {
      return LocalFunctionType( entity, discreteFunction_ );
    }

    DiscreteFunctionType &function ()
    {
      return discreteFunction_;
    }

  private:
    LocalFunctionStorage ( const ThisType & );
    ThisType operator= ( const ThisType & );

    DiscreteFunctionType &discreteFunction_;
  };



  namespace
  {

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

  template< class FunctionImp, class GridPartImp >
  class ConvertToGridFunction
  : public Fem::Function< typename FunctionImp::FunctionSpaceType,
                          ConvertToGridFunction< FunctionImp, GridPartImp > >,
    public HasLocalFunction
  {
    typedef ConvertToGridFunction< FunctionImp, GridPartImp > ThisType;
    typedef Fem::Function< typename FunctionImp::FunctionSpaceType, ThisType > BaseType;

    static const bool hasLocalFunction
      = Conversion< FunctionImp, HasLocalFunction >::exists;

    typedef ConvertDFTypeHelper< FunctionImp, GridPartImp, hasLocalFunction >
      Helper;
    typedef typename Helper::FunctionType ConvertedType;

  public:
    typedef FunctionImp FunctionType;
    typedef GridPartImp GridPartType;

    //! type of discrete function space 
    typedef typename ConvertedType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    // type of discrete function space
    typedef typename ConvertedType::FunctionSpaceType FunctionSpaceType; 

    //! type of grid 
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

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
    ConvertToGridFunction ( const std::string &name,
                            const FunctionImp &function,
                            const GridPartType &gridPart )
    : name_( name ),
      helper_( name, function, gridPart )
    {}

    ConvertToGridFunction ( const ThisType &other ) 
    : name_( other.name_ ),
      helper_( other.helper_ )
    {}

    //! evaluate function on local coordinate local 
    void evaluate ( const DomainType &global, RangeType &result ) const 
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
    const std::string &name() const
    {
      return name_;
    }

    const DiscreteFunctionSpaceType &space() const
    {
      return helper_.space();
    }

  private:
    const std::string name_;
    Helper helper_;
  };

  template< class Function, class GridPart >
  inline ConvertToGridFunction< Function, GridPart >
  convertToGridFunction ( const std::string &name,
                          const Function &function,
                          const GridPart &gridPart )
  {
    return ConvertToGridFunction< Function, GridPart >( name, function, gridPart );
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
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
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
  template< class EvalImp >
  class LocalFunctionAdapter
  : public Fem::Function< typename EvalImp::FunctionSpaceType, LocalFunctionAdapter< EvalImp > >,
    public HasLocalFunction
  {
    typedef LocalFunctionAdapter< EvalImp > ThisType;
    typedef Fem::Function< typename EvalImp::FunctionSpaceType, ThisType > BaseType;

  public:  
    //! Evaluate class
    typedef EvalImp EvalType;

    //! type of function 
    typedef typename BaseType::FunctionType FunctionType;

    //! traits class
    typedef LocalFunctionAdapterTraits< EvalType > Traits;

    //! type of grid part 
    typedef typename EvalImp::GridPartType GridPartType;
    
    //! type of discrete function
    typedef DiscreteFunctionSpaceAdapter< typename EvalType::FunctionSpaceType,
                                          typename EvalType::GridPartType >
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

    //! constructer taking instance of EvalImp class 
    LocalFunctionAdapter ( const std::string &name,
                           EvalType &eval,
                           const GridPartType &gridPart,
                           unsigned int order = DiscreteFunctionSpaceType::polynomialOrder )
    : space_( gridPart, order ),
      eval_( eval ),
      name_( name )
    {}

    // reference to function this local belongs to
    LocalFunctionAdapter( const ThisType &other )
    : space_( other.space_ ),
      eval_( other.eval_ ),
      name_( other.name_ )
    {}

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

  private:    
    DiscreteFunctionSpaceType space_; 
    mutable EvalType& eval_;
    const std::string name_;
  };


  template< class EvalImp >
  struct LocalFunctionAdapter< EvalImp >::LocalFunction
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

    //! constructor initializing local function 
    LocalFunction(const EntityType& en, const ThisType& a)
    : eval_(a.eval_)
    {
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
    void evaluate ( const PointType &x, RangeType &ret ) const
    {
      eval_.evaluate(x,ret);
    }

    //! jacobian of local function 
    template< class PointType >
    void jacobian ( const PointType &x, JacobianRangeType &ret ) const
    {
      eval_.jacobian( x, ret );
    }

    //! init local function
    void init(const EntityType& en)
    {
      eval_.init(en);
    } 
  private:
    EvalType& eval_;
  };

} // end namespace Dune

//@}

#endif // #ifndef DUNE_DISCRETEFUNCTIONADAPTER_HH
