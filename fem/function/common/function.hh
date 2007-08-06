#ifndef DUNE_FUNCTION_HH
#define DUNE_FUNCTION_HH

//- Dune includes 
#include <dune/common/fvector.hh>

//- local includes 
#include <dune/fem/operator/common/mapping.hh>

namespace Dune
{
/** @addtogroup Functions
    Functions are mapping from one finite dimensional
    vector space into another, e.g.,
    \f$K^n\f$ into \f$L^m\f$.
    They are element of a 
    \ref FunctionSpaceInterface "function space".

    \remark The interface is given by Function.
    @{
**/
  
/** \brief 
    Abstract class representing a function

    Template parameters are:
    -  FunctionSpaceImp      type of the function space where the function 
                             belongs to.
    -  FunctionImp           type of the implemented function (Barton-Nackman)

    \interfaceclass
*/
template< class FunctionSpaceImp, class FunctionImp>
class Function : 
    public Mapping < typename FunctionSpaceImp::DomainFieldType,
                     typename FunctionSpaceImp::RangeFieldType ,
                     typename FunctionSpaceImp::DomainType, 
                     typename FunctionSpaceImp::RangeType > 
{
  //! type of this 
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

  //! constructor
  inline explicit Function ( const FunctionSpaceType &fSpace )
  : functionSpace_( fSpace )
  {
  }   

  //! copy constructor
  Function (const FunctionType& org) : functionSpace_ (org.functionSpace_) 
  {
  }   

  //! destructor 
  virtual ~Function () {} 

  /** \brief application operator call evaluate 
      \param[in] arg argument 
      \param[out] dest destination, i.e. f(arg) 
  */
  virtual void operator()(const DomainType & arg, RangeType & dest) const 
  {
    evaluate(arg,dest);
  }

  /** \brief evaluate function f 
      \param[in] arg argument 
      \param[out] dest destination, i.e. f(arg) 
  */
  void evaluate(const DomainType & arg, RangeType & dest) const 
  {
    asImp().evaluate(arg, dest);
  }

  /** \brief evaluate function and derivatives specialization via derivation: 
      \param[in] diffVariable derivation 0 == evaluate function, derivation 1 == evaluate partial derivative diffVariable[0] (0 == x0, 1 == x1,..., n == xn)
      \param[in] arg argument 
      \param[out] dest destination, i.e. f(arg) 
  */    
  template <int derivation>
  void evaluate  ( const FieldVector<deriType, derivation> &diffVariable, 
                   const DomainType& arg, RangeType & dest) const { 
    asImp().evaluate(diffVariable, arg, dest);
  }

  /** \brief Get access to the related function space
      \return return reference to function space 
  */ 
  const FunctionSpaceType& space() const { return functionSpace_; }

protected:
  //! Barton-Nackman trick
  FunctionImp& asImp() { 
    return static_cast<FunctionImp&>(*this); 
  }
  //! Barton-Nackman trick
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
  virtual void apply (const DomainType& arg, RangeType& dest) const 
  {
    operator()(arg, dest);
  }
};

///@} 
}
#endif
