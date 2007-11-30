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
      i.e., expecting \ref LocalFunctionInterface "local functions"
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

      //! evaluate local function 
      template< class PointType >
      void evaluate ( const PointType &x, RangeType &ret ) const
      {
        DomainType global = geometry_->global( coordinate( x ) );
        function_.evaluate( global, ret );
      }

      //! evaluate local function 
      template< class QuadratureType >
      void evaluate( const QuadratureType &quadrature,
                     const int quadPoint,
                     RangeType &ret ) const 
      {
        evaluate( quadrature[ quadPoint ], ret );
      }

      //! jacobian of local function 
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        DomainType global = geometry_->global( coordinate( x ) );
        function_.jacobian( global, ret );
      }

      //! jacobian of local function 
      template< class QuadratureType >
      void jacobian ( const QuadratureType &quadrature,
                      const int quadPoint,
                      JacobianRangeType &ret ) const 
      {
        jacobian( quadrature[ quadPoint ], ret );
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

    // reference to function this local belongs to
    inline DiscreteFunctionAdapter
      ( const std :: string &name,
        const FunctionType &f,
        const GridPartType &gridPart,
        unsigned int order = DiscreteFunctionSpaceType :: polynomialOrder )
    : BaseType(space_),
      space_( gridPart, order ),
      function_( f ),
      name_( name )
    {
    }

    // reference to function this local belongs to
    DiscreteFunctionAdapter( const ThisType &other ) 
    : BaseType( other ),
      space_( other.space_ ),
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
    //! reference to function 
    const FunctionType& function_; 
    
    const std::string name_;
  };

} // end namespace Dune

//@}

#endif
