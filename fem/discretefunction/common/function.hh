#ifndef DUNE_FUNCTION_HH
#define DUNE_FUNCTION_HH

//- Dune includes 
#include <dune/common/fvector.hh>

//- local includes 
#include "../../operator/common/mapping.hh"

namespace Dune{
/** @defgroup FunctionCommon Functions
    @ingroup Functions
    Functions are Mappings from \f$K^n\f$ into \f$L^m\f$ where 
    \f$K\f$ and \f$L\f$ are fields.
    @{
 */

//! type of derivative identifier 
typedef int deriType;
  
/*!  Abstract class representing a function

    Template parameters are:
    -  FunctionSpaceImp      type of the function space where the function 
                             belongs to.
    -  FunctionImp           type of the implemented function (Barton-Nackman)
*/
template< class FunctionSpaceImp, class FunctionImp>
class Function : 
    public Mapping < typename FunctionSpaceImp::DomainFieldType,
                     typename FunctionSpaceImp::RangeFieldType ,
                     typename FunctionSpaceImp::DomainType, 
                     typename FunctionSpaceImp::RangeType > 
{

  // type of this 
  typedef Function<FunctionSpaceImp,FunctionImp> FunctionType;

public:
  //! type of function space this function belongs to 	
  typedef FunctionSpaceImp FunctionSpaceType;
  //! domain type (from function space)
  typedef typename FunctionSpaceType::DomainType DomainType ;
  //! range type (from function space)
  typedef typename FunctionSpaceType::RangeType RangeType ;
  //! jacobian type (from function space)
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  //! hessian type (from function space)
  typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;

  //! Constructor
  Function (const FunctionSpaceType & f) : functionSpace_ (f) {} ;  

  //! application operator
  virtual void operator()(const DomainType & arg, RangeType & dest) const 
  {
    evaluate(arg,dest);
  }

  //! evaluate Function
  void evaluate(const DomainType & arg, RangeType & dest) const 
  {
    asImp().evaluate(arg, dest);
  }

  //! evaluate Function
  void eval(const DomainType & arg, RangeType & dest) const DUNE_DEPRECATED  {
    asImp().eval(arg, dest);
  }

  /** \brief evaluate function and derivatives 
     specialization via derivation: 
	derivation 0 == evaluate function 
	derivation 1 == evaluate partial derivative diffVariable[0] (0 == x0, 1 == x1,..., n == xn)
	*/		
  template <int derivation>
  void evaluate  ( const FieldVector<deriType, derivation> &diffVariable, 
                   const DomainType& arg, RangeType & dest) const { 
    asImp().evaluate(diffVariable, arg, dest);
  }

  //! Get access to the related function space
  const FunctionSpaceType& getFunctionSpace() const DUNE_DEPRECATED { return functionSpace_; }

  //! Get access to the related function space
  const FunctionSpaceType& space() const { return functionSpace_; }

protected:
  //! Barton-Nackman trick
  FunctionImp& asImp() { 
    return static_cast<FunctionImp&>(*this); 
  }
  const FunctionImp& asImp() const { 
    return static_cast<const FunctionImp&>(*this); 
  }
  
  //! The related function space
  const FunctionSpaceType & functionSpace_;
private:
  //! Helper function for Mapping
  //! With this function, a combined mapping can choose the right application
  //! operator (i.e. the one from Mapping itself, or from Function/Operator)
  //! \note: Do not override this definition
  virtual void apply (const DomainType& arg, RangeType& dest) const {
    operator()(arg, dest);
  }
};


//! adapter to provide local function for a function
template <class FunctionImp, class GridType>
class LocalFunctionAdapter
{
  // type of function 
  typedef FunctionImp FunctionType;

  //! domain type (from function space)
  typedef typename FunctionType::DomainType DomainType ;
  //! range type (from function space)
  typedef typename FunctionType::RangeType RangeType ;
  //! jacobian type (from function space)
  typedef typename FunctionType::JacobianRangeType JacobianRangeType;

  //! type of codim 0 entity
  typedef typename GridType :: template Codim<0> :: Entity EntityType; 

  private:
  class LocalFunction
  {
    //! type of geometry 
    typedef typename EntityType :: Geometry GeometryImp;
  public:  
    //! constructor initializing local function 
    LocalFunction(const FunctionType& function, 
                  const EntityType& en)
      : function_(function) 
      , geometry_(&(en.geometry())) 
    {}

    //! copy constructor 
    LocalFunction(const LocalFunction& org) 
      : function_(org.function_) 
      , geometry_(org.geometry_)  
    {}

    //! evaluate local function 
    void evaluate(const DomainType& local, RangeType& result) const
    {
      DomainType global = geometry_->global(local);
      function_.evaluate(global,result);
    }

    //! evaluate local function 
    template <class QuadratureType>
    void evaluate(const QuadratureType& quad,
                  const int quadPoint, 
                  RangeType& result) const 
    {
      function_.evaluate(quad.point(quadPoint), result);
    }

    //! evaluate local function 
    void jacobian(const DomainType& local, RangeType& result) const
    {
      assert(false);
      DomainType global = geometry_->global(local);
      function_.evaluate(global,result);
    }

    //! evaluate local function 
    template <class QuadratureType>
    void jacobian(const QuadratureType& quad,
                  const int quadPoint, 
                  RangeType& result) const 
    {
      assert(false);
      function_.evaluate(quad.point(quadPoint), result);
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
  LocalFunctionAdapter(const FunctionType& f) 
    : function_(f)
  {}

  //! evaluate function on local coordinate local 
  void evaluate(const DomainType& global, RangeType& result) const 
  {
    function_.evaluate(global,result);  
  }

  //! return local function object 
  LocalFunctionType localFunction(const EntityType& en) const 
  {
    return LocalFunctionType(function_,en);
  }
private:    
  //! reference to function 
  const FunctionType& function_; 
};

/** @} end documentation group */

}
#endif
