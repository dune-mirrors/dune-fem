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
    typedef Function<DiscreteFunctionSpaceAdapter<typename FunctionImp :: FunctionSpaceType,GridPartImp>, 
                     DiscreteFunctionAdapter<FunctionImp,GridPartImp> > BaseType;
  public:  
    //! type of grid part 
    typedef GridPartImp GridPartType;
                       
    //! type of traits 
    typedef DiscreteFunctionAdapterTraits<FunctionImp,GridPartImp> Traits;
    //! type of this pointer 
    typedef DiscreteFunctionAdapter<FunctionImp,GridPartImp> ThisType;

    // type of function 
    typedef FunctionImp FunctionType;

    // type of discrete function space
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType; 

    //! type of discrete function space 
    typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

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
        , global_(0)
      {}

      LocalFunction(const ThisType& a)
        : function_(a.function_) 
        , geometry_(0) 
        , global_(0)
      {}

      //! copy constructor 
      LocalFunction(const LocalFunction& org) 
        : function_(org.function_) 
        , geometry_(org.geometry_)  
        , global_(0)
      {}

      //! evaluate local function 
      void evaluate(const DomainType& local, RangeType& result) const
      {
        global_ = geometry_->global(local);
        function_.evaluate(global_,result);
      }

      //! evaluate local function 
      template <class QuadratureType>
      void evaluate(const QuadratureType& quad,
                    const int quadPoint, 
                    RangeType& result) const 
      {
        evaluate(quad.point(quadPoint), result);
      }

      //! jacobian of local function 
      void jacobian(const DomainType& local, JacobianRangeType& result) const
      {
        assert(false);
        abort();
      }

      //! jacobian of local function 
      template <class QuadratureType>
      void jacobian(const QuadratureType& quad,
                    const int quadPoint, 
                    JacobianRangeType& result) const 
      {
        assert(false);
        abort();
      }

      //! init local function
      void init(const EntityType& en) 
      {
        geometry_ = &(en.geometry());
      } 

    private:
      const FunctionType& function_;
      const GeometryImp* geometry_;
      mutable DomainType global_;
    };

    public:
    //! type of local function to export 
    typedef LocalFunction LocalFunctionType; 

    // reference to function this local belongs to
    DiscreteFunctionAdapter(const std::string name, const FunctionType& f, const GridPartType& gridPart) 
      : BaseType(space_)
      , space_(gridPart)
      , function_(f)
      , name_(name)
    {}

    // reference to function this local belongs to
    DiscreteFunctionAdapter(const DiscreteFunctionAdapter& org) 
      : BaseType(org)
      , space_(org.space_)
      , function_(org.function_)
      , name_(org.name_)
    {}

    //! evaluate function on local coordinate local 
    /** \brief @copydoc DiscreteFunctionInterface::localFunction */ 
    void evaluate(const DomainType& global, RangeType& result) const 
    {
      function_.evaluate(global,result);  
    }

    /** \brief @copydoc DiscreteFunctionInterface::localFunction */ 
    const LocalFunctionType localFunction(const EntityType& en) const 
    {
      return LocalFunctionType(en,*this);
    }

    /** \brief @copydoc DiscreteFunctionInterface::localFunction */ 
    LocalFunctionType localFunction(const EntityType& en) 
    {
      return LocalFunctionType(en,*this);
    }

    /** \brief @copydoc DiscreteFunctionInterface::name */ 
    const std::string& name() const 
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
#include "discretefunction.cc"
#endif
