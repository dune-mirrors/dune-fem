#ifndef DUNE_FUNCTIONSPACEINTERFACE_HH
#define DUNE_FUNCTIONSPACEINTERFACE_HH

#include <dune/common/fvector.hh>

namespace Dune{

/*! @addtogroup FunctionSpace 
  This provides the interfaces for analytical function spaces. 
  A function space is characterized by it's domain and range field type and  
  the two finite dimensional vector spaces over these fields. 
  
  All corresponding types 
  are given in a Traits class. In addition to the
  domain and range types also types for 
  the jacobian and the hessian must be provieded.

  \remarks The interface for analytical functions spaces 
  is defined by the class FunctionSpaceInterface.
  @{
 */

/*! \brief An arbitrary function space
    Base class for specific function spaces.

    \interfaceclass
*/
template<typename FunctionSpaceTraits> 
class FunctionSpaceInterface {
public:
/*! Dimensions of the domain and range field */
  enum { DimDomain = FunctionSpaceTraits::DimDomain ,//!< dimension of range vector space 
         DimRange = FunctionSpaceTraits::DimRange    //!< dimension of domain vector space 
  };
  /** \brief Intrinsic type used for values in the domain field (usually a double) */
  typedef typename FunctionSpaceTraits::DomainFieldType DomainFieldType;
  
  /** \brief Intrinsic type used for values in the range field (usually a double) */
  typedef typename FunctionSpaceTraits::RangeFieldType RangeFieldType;
  
  /** \brief Type of domain vector (using type of domain field) 
      has a Dune::FieldVector type interface */
  typedef typename FunctionSpaceTraits::DomainType DomainType;
  
  /** \brief Type of range vector (using type of range field) 
      has a Dune::FieldVector type interface */
  typedef typename FunctionSpaceTraits::RangeType RangeType;
  
  /** \brief Intrinsic type used for the jacobian values 
      has a Dune::FieldMatrix type interface */
  typedef typename FunctionSpaceTraits::LinearMappingType JacobianRangeType;
  
  /** \brief Intrinsic type used for the hessian values 
      has a Dune::FieldMatrix type interface */
  typedef FieldVector<JacobianRangeType, DimRange> HessianRangeType;
  
  /** \brief Type of corresponding scalar space 
      is of FunctionSpaceInterface type with scalar 
      range */
  typedef typename FunctionSpaceTraits::ScalarFunctionSpaceType ScalarFunctionSpaceType;
};


///@}
} // end namespace Dune 
#endif
