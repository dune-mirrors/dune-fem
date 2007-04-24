#ifndef DUNE_LOCALFUNCTION_HH
#define DUNE_LOCALFUNCTION_HH

namespace Dune {


/** @defgroup LocalFunction LocalFunction 
   @ingroup DiscreteFunction

  On every element from a discrete function the local funtion can be accessed.
  With the local function one has access to the dof and on the other hand to 
  the base function set of this actual element. Therefore this is called a 
  local function. 

  @{
 */

  
//****************************************************************************
//
//  --LocalFunctionInterface 
//
//! The LocalFunctionInterface is the Interface to local function which
//! form the discrete Function 
//
//****************************************************************************
template < class DiscreteFunctionSpaceType, class LocalFunctionImp > 
class LocalFunctionInterface 
{
public:
  //! this are the types for the derived classes 
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  //! type of base function set  
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType; 

  //! access to dof number num, all dofs of the local function
  RangeFieldType& operator [] (int num) 
  {
    return asImp().operator [] ( num );
  }

  //! access to dof number num, all dofs of the local function
  const RangeFieldType& operator [] (int num) const 
  {
    return asImp().operator [] ( num );
  }

  //! return the number of local dof of this local function 
  int numDofs() const 
  {
    return asImp().numDofs();
  }
  
  //! evaluate local function.
  //! gets an x in local coordinates
  void evaluate (const DomainType & x, RangeType & ret)
  {
    asImp().evaluate(x,ret);
  }

  //! evaluate jacobian on reference element coordinate x
  template <class EntityType>
  void jacobian(EntityType& en, 
		const DomainType& x, 
		JacobianRangeType& ret) 
  {
    asImp().jacobian(en,x,ret);
  }

  //! return reference to corresponding base function set of local
  //! function
  const BaseFunctionSetType& baseFunctionSet() const 
  {
    return asImp().baseFunctionSet();
  }
   
  void assign(int dofNum, const RangeType& dofs) DUNE_DEPRECATED {
    asImp().assign(dofNum, dofs);
  }
private:
  //! Barton-Nackman trick 
  LocalFunctionImp& asImp() 
  { 
    return static_cast<LocalFunctionImp&>(*this);
  }
  
  const LocalFunctionImp& asImp() const  
  {
    return static_cast<const LocalFunctionImp&>(*this);
  }
}; // end LocalFunctionInterface

//************************************************************************
//
//  --LocalFunctionDefault 
//
//! The Interface to the dune programmer, use this class to derive 
//! the own implementation. But only the methods declared in the interface
//! class must be implemented. 
//!
//************************************************************************
template < class DiscreteFunctionSpaceType, class LocalFunctionImp > 
class LocalFunctionDefault : public LocalFunctionInterface <
DiscreteFunctionSpaceType , LocalFunctionImp > 
{
public: 
  //! these are the types for the derived classes 
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
 
public:
  //! Constructor
  LocalFunctionDefault() : xLoc_(0.0) {}

  //! evaluate the local function on real world coordinate x and return ret 
  template <class EntityType>
  void evaluateGlobal(EntityType& en, 
                      const DomainType& x, 
                      RangeType& ret)
  {
    ret = 0.0;
    xLoc_ = en.geometry().local(x);
    evaluateLocal(en,xLoc_,ret);
  }

  //! Evaluation using a quadrature
  template <class QuadratureType>
  void evaluate(QuadratureType& quad,
                int quadPoint,
                RangeType& ret)
  {
    evaluate(quad.point(quadPoint), ret);
  }
  
  //! jacobian of the local function using real world coordinate x
  template <class EntityType>
  void jacobianGlobal(EntityType& en, 
                      const DomainType& x, 
                      JacobianRangeType& ret) 
  {
    ret = 0.0;
    xLoc_ = en.geometry().local(x);
    jacobianLocal(en,xLoc_,ret);
  }

  //! Evaluation of jacobian using a quadrature
  template <class EntityType, class QuadratureType>
  void jacobian(EntityType& en,
                QuadratureType& quad,
                int quadPoint,
                JacobianRangeType& ret)
  {
    jacobian(en, quad.point(quadPoint), ret);
  }

private:
  mutable DomainType xLoc_;

}; // end LocalFunctionDefault

} // end namespace Dune 

#endif

