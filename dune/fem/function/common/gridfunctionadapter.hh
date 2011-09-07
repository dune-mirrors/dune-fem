#ifndef DUNE_FEM_GRIDFUNCTIONADAPTER_HH
#define DUNE_FEM_GRIDFUNCTIONADAPTER_HH

//- local includes 
#include <dune/fem/version.hh>
#include <dune/fem/function/common/discretefunction.hh>

// for compatibility 
#include <dune/fem/function/common/localfunctionadapter.hh>

namespace Dune
{

  /** @addtogroup GridFunctionAdapter

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
  class GridFunctionAdapter;

  //! traits of GridFunctionAdapter 
  template <class FunctionImp, class GridPartImp> 
  struct GridFunctionAdapterTraits 
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

    typedef GridFunctionAdapter<FunctionImp,GridPartImp> DiscreteFunctionType;
  };



  // GridFunctionAdapter
  // -----------------------

  /** \brief GridFunctionAdapter provides local functions for a Function. 
   */
  template< class FunctionImp, class GridPartImp >
  class GridFunctionAdapter 
  : public Fem::Function< typename FunctionImp::FunctionSpaceType,
                          GridFunctionAdapter< FunctionImp, GridPartImp > >,
    public HasLocalFunction
  {
    typedef GridFunctionAdapter< FunctionImp, GridPartImp > ThisType;
    typedef Fem::Function< typename FunctionImp::FunctionSpaceType, ThisType > BaseType;

    // Make sure the function is not a discrete functon
    dune_static_assert( !(Conversion< FunctionImp, HasLocalFunction >::exists),
                        "FunctionType may not be a discrete function type." );

  public:  
    //! type of traits 
    typedef GridFunctionAdapterTraits< FunctionImp, GridPartImp > Traits;

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
    GridFunctionAdapter ( const std::string &name,
                          const FunctionType &f,
                          const GridPartType &gridPart,
                          unsigned int order = DiscreteFunctionSpaceType::polynomialOrder )
    : space_( gridPart, order ),
      localFunctionStorage_( *this ),
      function_( f ),
      name_( name )
    {}

    // reference to function this local belongs to
    GridFunctionAdapter( const ThisType &other ) 
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
    //! evaluate function on local coordinate local 
    void jacobian(const DomainType& global, JacobianRangeType& result) const 
    {
      function_.jacobian(global,result);  
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



  // GridFunctionAdapter::LocalFunction
  // --------------------------------------

  template< class Function, class GridPart >
  class GridFunctionAdapter< Function, GridPart >::LocalFunction
  {
    typedef LocalFunction ThisType;
    typedef GridFunctionAdapter< Function, GridPart > DiscreteFunctionType;

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



  // GridFunctionAdapter::LocalFunctionStorage
  // ---------------------------------------------

  template< class Function, class GridPart >
  class GridFunctionAdapter< Function, GridPart >::LocalFunctionStorage
  {
    typedef LocalFunctionStorage ThisType;
    typedef GridFunctionAdapter< Function, GridPart > DiscreteFunctionType;

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


  //- deprecated class 
  template< class FunctionImp, class GridPartImp >
  class DiscreteFunctionAdapter 
    : public GridFunctionAdapter< FunctionImp, GridPartImp > 
  {
    typedef GridFunctionAdapter< FunctionImp, GridPartImp > BaseType;
  public:  
    DUNE_VERSION_DEPRECATED(1,3,GridFunctionAdapter)
    DiscreteFunctionAdapter ( const std::string &name,
                              const FunctionImp &f,
                              const GridPartImp &gridPart,
                              unsigned int order = BaseType::DiscreteFunctionSpaceType::polynomialOrder )
    : BaseType( name, f, gridPart, order )
    {}
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
      GridFunctionAdapter<FunctionImp,GridPartType>
    {
      typedef ConvertDFTypeHelper<FunctionImp,GridPartType,false>
        ThisType;
      typedef GridFunctionAdapter<FunctionImp,GridPartType> BaseType;
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

} // end namespace Dune

//@}

#endif // #ifndef DUNE_DISCRETEFUNCTIONADAPTER_HH