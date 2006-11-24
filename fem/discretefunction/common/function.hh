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
                     typename FunctionSpaceImp::RangeType > {

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

/** @} end documentation group */

}

#endif
