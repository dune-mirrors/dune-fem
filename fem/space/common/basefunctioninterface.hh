#ifndef DUNE_BASEFUNCTIONINTERFACE_HH
#define DUNE_BASEFUNCTIONINTERFACE_HH

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>

//- local includes 
#include <dune/fem/operator/common/mapping.hh>

namespace Dune{

/** @addtogroup BaseFunction  
    
  The base functions are essential to describe a numerical solutions.
  Here the interface of base functions and the corresponding base 
  function set is presented. 
  The user always works with the
  BaseFunctionSets, where all different base functions for 
  on element type are known.

  \remarks  
  The interface for the BaseFunctionSet is defined by
  BaseFunctionSetInterface. 

  \remarks 
  The interface for using a BaseFunction is described by BaseFunctionInterface.
  Note that the member methods are virtual.
  @{
 */
  


/** \brief interface for base functions
 *
 *  A base function can be evaluated on a point from the Domain and the
 *  outcome is a point from Range. The Types of Domain and Range are stored
 *  by typedefs in FunctionSpaceType which is the template parameter of
 *  BaseFunctionInterface. 
 *  
 *  \remark The BaseFunctionInterface is (unlike most other interfaces in
 *          dune-fem) a virtual interface, i.e., all methods are declared
 *          virtual. The Barton-Nackman trick is not used. Therefore, all base
 *          base functions of a discrete function space are derived from the
 *          same type (all template parameters are the same).
 *
 *  \remark The virtual calls are usually not a problem, because the value of
 *          base functions in quadrature points is cached when using the
 *          CachingStorage.
 */
template< class FunctionSpaceImp >
class BaseFunctionInterface
: public Mapping< typename FunctionSpaceImp :: DomainFieldType,
                  typename FunctionSpaceImp :: RangeFieldType, 
                  typename FunctionSpaceImp :: DomainType, 
                  typename FunctionSpaceImp :: RangeType > 
{
    
public:
  //! \brief type of function space 
  typedef FunctionSpaceImp FunctionSpaceType;
  //! \brief type of domain vector, i.e. element of \f$R^{dimDomain}\f$
  typedef typename FunctionSpaceType::DomainType DomainType;
  //! \brief type of range vector, i.e. element of \f$R^{dimRange}\f$
  typedef typename FunctionSpaceType::RangeType RangeType;
  //! \brief type of domain vector components, i.e. double 
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  //! \brief type of range vector components, i.e. double 
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  //! \brief dimension of domain 
  enum { dimDomain = FunctionSpaceType::dimDomain };
  //! \brief dimension of range 
  enum { dimRange  = FunctionSpaceType::dimRange };

  //! empty constructor 
  BaseFunctionInterface () {}

  //! empty destructor 
  virtual ~BaseFunctionInterface() {}

  /** \brief 
     evaluate the function at Domain x, and store the value in Range Phi
     diffVariable stores information about which gradient is to be
     evaluated. In the derived classes the evaluate method are template
     methods with template parameter "length of Vec". 
     Though the evaluate Methods can be spezialized for each
     differentiation order
     \param x The local coordinate in the reference element
     \param diffVariable description of the derivative to be returned
     \param phi return value
  */
  virtual void evaluate ( const FieldVector<deriType, 0> &diffVariable, 
                          const DomainType &x , RangeType &phi) const = 0;

  /** diffVariable contain the component of the gradient which is delivered.
      for example gradient of the basefunction x component ==>
      diffVariable[0] == 0, y component ==> diffVariable[0] == 1 ...
  */
  virtual void evaluate ( const FieldVector<deriType, 1> &diffVariable, 
                          const DomainType & , RangeType &) const = 0; 

  /** \brief diffVariable contain the component of the hessian which is delivered.
      for example hessian of the basefunction x component ==>
      diffVariable[0][0] == 0, y component ==> diffVariable[0][1] == 1 ...
  */
  virtual void evaluate ( const FieldVector<deriType, 2> &diffVariable, 
                          const DomainType & , RangeType &) const = 0;

  //! \todo who has this inserted? 
  virtual int order() const { return -1; }
};

///@}
} // end namespace Dune 
#include "basefunctionsetinterface.hh"
#endif
