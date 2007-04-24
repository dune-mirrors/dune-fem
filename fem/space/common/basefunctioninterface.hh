#ifndef DUNE_BASEFUNCTIONINTERFACE_HH
#define DUNE_BASEFUNCTIONINTERFACE_HH

//- Dune includes 

//- local includes 
#include <dune/fem/operator/common/mapping.hh>

namespace Dune{

/** @defgroup BaseFunction Base Functions
  @ingroup DiscreteFunctionSpace
    
  The base functions are essential to describe a numerical solutions.
  Here the interface of base functions and the corresponding base 
  function set is presented. The user always works with the base function
  set, where all diffrent base functions for on element type are known.

  @{
 */

// just to make it easy to change 
typedef int deriType;

// just for make changes easy 
template <int dim>
struct DiffVariable
{
    typedef FieldVector<deriType, dim> Type;
};

//*************************************************************************
//! BaseFunctionInterface is the interface to a base function. 
//! A base function can be evaluated on a point from the Domain and the
//! outcome is a point from Range. The Types of Domain and Range are stored
//! by typedefs in FunctionSpaceType which is the template parameter of
//! BaseFunctionInterface. 
//*************************************************************************
template<class FunctionSpaceImp>
class BaseFunctionInterface 
: public Mapping< typename FunctionSpaceImp::DomainFieldType,
                  typename FunctionSpaceImp::RangeFieldType, 
                  typename FunctionSpaceImp::DomainType, 
                  typename FunctionSpaceImp::RangeType > 
{
    
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  enum { DimDomain = FunctionSpaceType::DimDomain };
  enum { DimRange  = FunctionSpaceType::DimRange  };

  BaseFunctionInterface () {}

  virtual ~BaseFunctionInterface() {}

  //! evaluate the function at Domain x, and store the value in Range Phi
  //! diffVariable stores information about which gradient is to be
  //! evaluated. In the derived classes the evaluate method are template
  //! methods with template parameter "length of Vec". 
  //! Though the evaluate Methods can be spezialized for each
  //! differentiation order
  //! \param x The local coordinate in the reference element
  //! \param diffVariable description of the derivative to be returned
  //! \param phi return value
  virtual void evaluate ( const FieldVector<deriType, 0> &diffVariable, 
                          const DomainType &x , RangeType &phi) const = 0;

  //! diffVariable contain the component of the gradient which is delivered.
  //! for example gradient of the basefunction x component ==>
  //! diffVariable[0] == 0, y component ==> diffVariable[0] == 1 ...
  virtual void evaluate ( const FieldVector<deriType, 1> &diffVariable, 
                          const DomainType & , RangeType &) const = 0; 

  //! diffVariable contain the component of the hessian which is delivered.
  //! for example hessian of the basefunction x component ==>
  //! diffVariable[0][0] == 0, y component ==> diffVariable[0][1] == 1 ...
  virtual void evaluate ( const FieldVector<deriType, 2> &diffVariable, 
                          const DomainType & , RangeType &) const = 0;

  virtual int order() const { return -1; }
};

/** @} end documentation group */


// @defgroup BaseFunctionSet BaseFunctionSet
/**  @ingroup BaseFunction
    
  The base functions are essential to describe a numerical solutions.
  Here the interface of base functions and the corresponding base 
  function set is presented. The user always works with the base function
  set, where all diffrent base functions for on element type are known.

//   @{
 */

//****************************************************************************
//
//  --BaseFunctionSetInterface
//
//! Why has the BaseFunctionInterface class virtual methods?
//!
//! Because the we want to have different base functions in 
//! at the same time and we havent found a way to do this with 
//! Barton-Nackman. But there is a solution which dosent cost to much
//! efficiency, for example in FastBaseFunctionSet all evaluations of the
//! BaseFunction are cached in a vector, though the virtual method evaluate
//! of the BaseFunctionImp has only to be called once for the same quadrature
//! rule. If you change the quadrature rule then on the first call the
//! evaluations are stored again.
//! This method brings us flexebility and effeciency. 
//!
//****************************************************************************
template<class BaseFunctionSetTraits> 
class BaseFunctionSetInterface  
{
public:
  typedef typename BaseFunctionSetTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename BaseFunctionSetTraits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;

  enum { DimDomain = FunctionSpaceType::DimDomain };
  enum { DimRange  = FunctionSpaceType::DimRange  };
  
  typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;

public:
  
  //! \todo Please doc me!
  BaseFunctionSetInterface () {}

  virtual ~BaseFunctionSetInterface() {}
  
  //! Number of base functions
  int numBaseFunctions() const {
    // stay on the safe side and call the deprecated version for now
    return asImp().numBaseFunctions();
  }

  //! \brief evaluate basefunction baseFunct for component diffVariable 
  //! on given point x, phi is the evaluation  
  template <int diffOrd>
  void evaluate (int baseFunct, 
                 const FieldVector<deriType, diffOrd> &diffVariable, 
                 const DomainType & x, RangeType & phi ) const 
  {
    asImp().evaluate( baseFunct, diffVariable, x, phi );
  }

  //! \brief evaluate basefunction baseFunct for component diffVariable 
  //! on quadrature point quadPoint, phi is the evaluation  
  template <int diffOrd, class QuadratureImp>
  void evaluate (int baseFunct, 
                 const FieldVector<deriType, diffOrd> &diffVariable, 
                 QuadratureImp& quad, 
                 int quadPoint, RangeType & phi ) const 
  {
    asImp().evaluate( baseFunct, diffVariable, quad, quadPoint, phi );
  }

private:
  //! Barton-Nackman trick 
  BaseFunctionSetType &asImp() { 
    return static_cast<BaseFunctionSetType&>(*this); 
  }
  
  const BaseFunctionSetType &asImp() const { 
    return static_cast<const BaseFunctionSetType&>(*this); 
  }

};

//*************************************************************************
//  
//  --BaseFunctionSetDefault
//
//! The BaseFunctionSetDefault class is the internal interface. Here some
//! default behavoir is implemented which always can be overloaded by the
//! implementation class, but not has to. 
//!
//*************************************************************************
template<class BaseFunctionSetTraits> 
class BaseFunctionSetDefault  : 
  public BaseFunctionSetInterface <BaseFunctionSetTraits> 
{
public:
  typedef typename BaseFunctionSetTraits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename BaseFunctionSetTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  enum { dimRow = JacobianRangeType::rows };
  enum { dimCol = JacobianRangeType::cols };
  
  typedef typename FunctionSpaceType::DomainType DomainType ;
  typedef typename FunctionSpaceType::RangeType RangeType ;
  typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;
  typedef typename FunctionSpaceType::RangeFieldType DofType;

public:
  //! set the default diffVar Types 
  BaseFunctionSetDefault () : 
    BaseFunctionSetInterface<BaseFunctionSetTraits> ()
  {
    for(int i=0; i<dimCol; i++)
      jacobianDiffVar_[i] = i;
  };
 
  virtual ~BaseFunctionSetDefault() {}

  //! default evaluate using the evaluate interface 
  void eval(int baseFunct, const DomainType & x, RangeType & phi) const DUNE_DEPRECATED
  {
    asImp().evaluate(baseFunct, diffVariable_ , x , phi);
    return;
  }

  //! default evaluate using the evaluate interface 
  void evaluate(int baseFunct, const DomainType & x, RangeType & phi) const
  {
    asImp().evaluate(baseFunct, diffVariable_ , x , phi);
    return;
  }

  //! default implementation for evaluation 
  template <class QuadratureImp>
  void eval(int baseFunct, QuadratureImp & quad, int quadPoint, RangeType & phi) const DUNE_DEPRECATED;

  //! default implementation for evaluation 
  template <class QuadratureImp>
  void evaluate(int baseFunct, QuadratureImp & quad, int quadPoint, RangeType & phi) const 
  {
    asImp().evaluate( baseFunct, diffVariable_ , quad, quadPoint, phi );
    return;
  }

  //! default evaluate using the evaluate interface 
  void jacobian(int baseFunct, const DomainType & x, JacobianRangeType & phi) const 
  {
    RangeType tmp;
    for(int i=0; i<dimCol; i++)
    {
      asImp().evaluate( baseFunct, jacobianDiffVar_[i] , x , tmp );
      for(int j=0; j<dimRow; j++)
        phi[j][i] = tmp[j];
    }
  }

  //! default implementation of evaluation the gradient 
  template <class QuadratureImp>
  void jacobian ( int baseFunct, QuadratureImp & quad, 
                  int quadPoint, JacobianRangeType & phi ) const 
  {
    RangeType tmp;
  
    for(int i=0; i<dimCol; i++)
    {
      asImp().evaluate( baseFunct, jacobianDiffVar_[i] , quad, quadPoint, tmp );
      for(int j=0; j<dimRow; j++)
        phi[j][i] = tmp[j];
    }
  }

  //! \brief evaluate basefunction and multiply with factor, return is DofType 
  DofType evaluateSingle(int baseFunct, 
                         const DomainType& xLocal,
                         const RangeType& factor) const
  {
    RangeType phi(0.);
    evaluate(baseFunct, xLocal, phi);
    return phi*factor;
  }

  //! \brief evaluate basefunction on quadrature point 
  //! and multiply with factor, return is DofType 
  template <class QuadratureType>
  DofType evaluateSingle(int baseFunct,
                         const QuadratureType& quad, int quadPoint,
                         const RangeType& factor) const
  {
    return evaluateSingle(baseFunct, quad.point(quadPoint), factor);
  }

  //! \brief evaluate gradient of basefunction on given entity (uses
  //! jacobianInverseTransposed) and multiply with factor, return is DofType 
  template <class Entity>
  DofType evaluateGradientSingle(int baseFunct,
                                 Entity& en,
                                 const DomainType& xLocal,
                                 const JacobianRangeType& factor) const
  {
    JacobianRangeType gradPhi(0.);
    jacobian(baseFunct, xLocal, gradPhi);

    DofType result = 0;
    for (int i = 0; i < FunctionSpaceType::DimRange; ++i) 
    {
      DomainType gradScaled(0.);
      en.geometry().jacobianInverseTransposed(xLocal).
        umv(gradPhi[i], gradScaled);
      result += gradScaled*factor[i];
    }
    return result;
  }

  //! \brief evaluate gradient of basefunction on quadrature point 
  //! on given entity (uses
  //! jacobianInverseTransposed) and multiply with factor, return is DofType 
  template <class Entity, class QuadratureType>
  DofType evaluateGradientSingle(int baseFunct,
                                 Entity& en,
                                 const QuadratureType& quad, int quadPoint,
                                 const JacobianRangeType& factor) const
  {
    return evaluateGradientSingle(baseFunct, en, quad.point(quadPoint),factor);
  }

private: 
  //! just diffVariable for evaluation of the functions 
  const FieldVector<deriType, 0> diffVariable_;

  //! diff variable for jacobian evaluation 
  FieldVector<deriType, 1> jacobianDiffVar_[dimCol];

  //! Barton-Nackman trick 
  BaseFunctionSetType &asImp() { return static_cast<BaseFunctionSetType&>(*this); }
  const BaseFunctionSetType &asImp() const 
  { return static_cast<const BaseFunctionSetType&>(*this); }

};

template<class BaseFunctionSetTraits> 
template <class QuadratureImp>
inline void BaseFunctionSetDefault<BaseFunctionSetTraits>::
eval(int baseFunct, QuadratureImp & quad, int quadPoint, RangeType & phi) const
{
  // call other eval because deprecated does not work for template methods 
  asImp().eval( baseFunct, quad.point(quadPoint) , phi );
  return;
}

/** @} end documentation group */

} // end namespace Dune 
#endif
