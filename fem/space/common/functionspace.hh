#ifndef DUNE_FUNCTIONSPACE_HH
#define DUNE_FUNCTIONSPACE_HH

#include <dune/common/fmatrix.hh>

namespace Dune{

/** @defgroup FunctionSpace FunctionSpace
  @ingroup Function
  This provides the interfaces for continuous function spaces. 
  A function space is characterized by it's domain and range field type and  
  the two vector spaces over these fields. All corresponding types 
  are given in a Traits class. 
  @{
 */

/** \brief An arbitrary function space
    Base class for specific function spaces.
*/
template<typename FunctionSpaceTraits> 
class FunctionSpaceBase {
public:
/** Dimensions of the domain and range field */
  enum { DimDomain = FunctionSpaceTraits::DimDomain, 
	 DimRange = FunctionSpaceTraits::DimRange};
/** Intrinsic type used for values in the domain field (usually a double) */
  typedef typename FunctionSpaceTraits::DomainFieldType DomainFieldType;
/** Intrinsic type used for values in the range field (usually a double) */
  typedef typename FunctionSpaceTraits::RangeFieldType RangeFieldType;
/** Type of domain vector (using type of domain field) */
  typedef typename FunctionSpaceTraits::DomainType DomainType;
/** Type of range vector (using type of range field) */
  typedef typename FunctionSpaceTraits::RangeType RangeType;
/** Intrinsic type used for the jacobian values */
  typedef typename FunctionSpaceTraits::LinearMappingType JacobianRangeType;
/** Intrinsic type used for the hessian values */
  typedef FieldVector<JacobianRangeType, DimRange> HessianRangeType;
/** Type of corresponding scalar space */
  typedef typename FunctionSpaceTraits::ScalarFunctionSpaceType ScalarFunctionSpaceType;
};
/** \brief Base class for vector valued function spaces.
*/
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
class FunctionSpace;
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
struct VectorSpaceTraits {
  typedef DomainFieldImp DomainFieldType;
  typedef RangeFieldImp RangeFieldType;
  typedef FieldVector<DomainFieldImp, n> DomainType;
  typedef FieldVector<RangeFieldImp, m> RangeType;
  typedef FieldMatrix<RangeFieldImp, m, n> LinearMappingType;
  typedef FunctionSpace<DomainFieldImp,RangeFieldImp,n,1> ScalarFunctionSpaceType;
  enum { DimDomain = n, DimRange = m};
};
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
class FunctionSpace : public
      FunctionSpaceBase<VectorSpaceTraits<DomainFieldImp,RangeFieldImp,n,m> > {};

/** \brief Base class for matrix valued function spaces.
*/
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
class MatrixFunctionSpace;

/** \brief RangeType class for matrix valued functions -
    derived from FieldMatrix but has representation as vector
*/
template <typename K,int n,int m> 
class RangeMatrix : public FieldMatrix<K,n,m> {
 public:
  typedef FieldMatrix<K,n,m> BaseType;
  enum {
    //! The number of rows.
    rows = BaseType::rows,
    //! The number of columns.
    cols = BaseType::cols,
    //! The total dimension of the matrix space.
    dimension = BaseType::rows*BaseType::cols
  };
  //===== constructors
  /** \brief Default constructor
   */
  RangeMatrix () : BaseType() {}
  
  /** \brief Constructor initializing the whole matrix with a scalar
   */
  RangeMatrix (const K& k) : BaseType(k) {}

  K& operator()(int r,int c) {
    return static_cast<BaseType&>(*this)[r][c];
  }
  const K operator()(int r,int c) const {
    return static_cast<BaseType&>(*this)[r][c];
  }
  K& operator[](int i) {
    int r = i/cols;
    int c = i%cols;
    return (*this)(r,c);
  }
  const K operator[](int i) const {
    int r = i/cols;
    int c = i%cols;
    return (*this)(r,c);
  }
  //! scalar product
  K operator* (const BaseType& y)
  {
    K ret(0);
    for (int i=0; i<n; i++)
      ret += static_cast<BaseType&>(*this)[i] * y[i];
    return ret;
  }

};
/** \brief JacobianRangeType class for matrix valued functions -
    derived from FieldMatrix 
*/
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
class MatrixMapping : 
    public FieldMatrix<RangeFieldImp,m1*m2,n> {
  // Implement L : VD -> VR where VR is the space of FieldMatrix<RFI,m1,m2>
  // VD is space of FieldVector<DFT,n>
 public:
  typedef FieldMatrix<RangeFieldImp,m1*m2,n> BaseType;
  typedef DomainFieldImp DomainFieldType;
  typedef RangeFieldImp RangeFieldType;
  typedef RangeMatrix<RangeFieldImp, m1,m2> RangeType;
  //===== constructors
  /** \brief Default constructor
   */
  MatrixMapping () : BaseType() {}
  /** \brief Constructor initializing the whole matrix with a scalar
   */
  MatrixMapping (const RangeFieldImp& k) : BaseType(k) {}
  
  FieldVector<DomainFieldImp,n>& operator[](int i) {
    return static_cast<BaseType&>(*this)[i];
  }
  const FieldVector<DomainFieldImp,n>& operator[](int i) const {
    return static_cast<BaseType&>(*this)[i];
  }
};
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
struct MatrixSpaceTraits {
  typedef DomainFieldImp DomainFieldType;
  typedef RangeFieldImp RangeFieldType;
  typedef FieldVector<DomainFieldImp, n> DomainType;
  typedef RangeMatrix<RangeFieldImp, m1,m2> RangeType;
  typedef MatrixMapping<DomainFieldImp,RangeFieldImp, n, m1,m2> LinearMappingType;
  typedef MatrixFunctionSpace<DomainFieldImp,RangeFieldImp,n,1,1> ScalarFunctionSpaceType;
  enum { DimDomain = n, DimRange = m1*m2};
};
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
class MatrixFunctionSpace : public
      FunctionSpaceBase<MatrixSpaceTraits<DomainFieldImp,RangeFieldImp,n,m1,m2> > {};

/** @} end documentation group */

}

#endif
