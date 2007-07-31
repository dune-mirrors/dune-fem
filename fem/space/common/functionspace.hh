#ifndef DUNE_FUNCTIONSPACE_HH
#define DUNE_FUNCTIONSPACE_HH

#include <dune/common/fmatrix.hh>

namespace Dune{

/*! @defgroup FunctionSpace FunctionSpace
  @ingroup FunctionCommon
  This provides the interfaces for continuous function spaces. 
  A function space is characterized by it's domain and range field type and  
  the two vector spaces over these fields. All corresponding types 
  are given in a Traits class. 
  @{
 */

/*! \brief An arbitrary function space
    Base class for specific function spaces.
*/
template<typename FunctionSpaceTraits> 
class FunctionSpaceBase {
public:
/*! Dimensions of the domain and range field */
  enum { DimDomain = FunctionSpaceTraits::DimDomain ,//!< dimension of range vector space 
         DimRange = FunctionSpaceTraits::DimRange    //!< dimension of domain vector space 
  };
  /** \brief Intrinsic type used for values in the domain field (usually a double) */
  typedef typename FunctionSpaceTraits::DomainFieldType DomainFieldType;
  
  /** \brief Intrinsic type used for values in the range field (usually a double) */
  typedef typename FunctionSpaceTraits::RangeFieldType RangeFieldType;
  
  /** \brief Type of domain vector (using type of domain field) */
  typedef typename FunctionSpaceTraits::DomainType DomainType;
  
  /** \brief Type of range vector (using type of range field) */
  typedef typename FunctionSpaceTraits::RangeType RangeType;
  
  /** \brief Intrinsic type used for the jacobian values */
  typedef typename FunctionSpaceTraits::LinearMappingType JacobianRangeType;
  
  /** \brief Intrinsic type used for the hessian values */
  typedef FieldVector<JacobianRangeType, DimRange> HessianRangeType;
  
  /** \brief Type of corresponding scalar space */
  typedef typename FunctionSpaceTraits::ScalarFunctionSpaceType ScalarFunctionSpaceType;
};

//! \brief Base class for vector valued function spaces.
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
class FunctionSpace;

//! \brief Traits class for vector function spaces 
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
struct VectorSpaceTraits {
  //! \brief @copydoc FunctionSpaceBase::DomainFieldType 
  typedef DomainFieldImp DomainFieldType;
  //! \brief @copydoc FunctionSpaceBase::RangeFieldType 
  typedef RangeFieldImp RangeFieldType;
  //! \brief @copydoc FunctionSpaceBase::DomainType 
  typedef FieldVector<DomainFieldImp, n> DomainType;
  //! \brief @copydoc FunctionSpaceBase::RangeType 
  typedef FieldVector<RangeFieldImp, m> RangeType;
  //! \brief linear mapping type  
  typedef FieldMatrix<RangeFieldImp, m, n> LinearMappingType;
  //! \brief scalar function space type 
  typedef FunctionSpace<DomainFieldImp,RangeFieldImp,n,1> ScalarFunctionSpaceType;
  //! \brief dimension of domain vector space 
  enum { DimRange = m};
  //! \brief dimension of range vector space 
  enum { DimDomain = n };
};

/** \brief FunctionSpace defines what the types of the domain vector
    space and the range vector space for a function are. 
*/
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
class FunctionSpace : public
      FunctionSpaceBase<VectorSpaceTraits<DomainFieldImp,RangeFieldImp,n,m> > {};

/*! \brief Base class for matrix valued function spaces.
*/
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
class MatrixFunctionSpace;

/*! \brief RangeType class for matrix valued functions -
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
  /*! \brief Default constructor
   */
  RangeMatrix () : BaseType() {}
  
  /*! \brief Constructor initializing the whole matrix with a scalar
   */
  RangeMatrix (const K& k) : BaseType(k) {}

  /** \brief access element in row r and column c 
      \param[in] r row 
      \param[in] c column
      \return reference to element in row r and column c 
  */
  K& operator()(int r,int c) {
    return static_cast<BaseType&>(*this)[r][c];
  }
  /** \brief access element in row r and column c 
      \param[in] r row 
      \param[in] c column
      \return reference to element in row r and column c 
  */
  const K operator()(int r,int c) const {
    return static_cast<BaseType&>(*this)[r][c];
  }
  
  /** \brief access i element where row = i/col and column = i%col 
      \param[in] i element number ot access 
      \return reference to element in row i/col and column i%col 
  */
  K& operator[](int i) {
    int r = i/cols;
    int c = i%cols;
    return (*this)(r,c);
  }
  /** \brief access i element where row = i/col and column = i%col 
      \param[in] i element number ot access 
      \return reference to element in row i/col and column i%col 
  */
  const K operator[](int i) const {
    int r = i/cols;
    int c = i%cols;
    return (*this)(r,c);
  }
  
  /** \brief scalar product
      \param y RangeMatrix to scalar multiply with 
      \return K scalar product 
  */
  K operator* (const BaseType& y)
  {
    K ret(0);
    for (int i=0; i<n; i++)
      ret += static_cast<BaseType&>(*this)[i] * y[i];
    return ret;
  }
  
  /** \brief vector space axpy operation
      \param a scalar factor 
      \param y RangeMatrix to multiply with  
      \return reference to this is returned (i.e. *this) 
  */
  RangeMatrix& axpy (const K& a, const BaseType& y)
  {
     for (int i=0; i<n; i++)
        static_cast<BaseType&>(*this)[i].axpy(a,y[i]);
     return *this;
  }
};
/*! \brief JacobianRangeType class for matrix valued functions -
    derived from FieldMatrix 
*/
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
class MatrixMapping : 
    public FieldMatrix<RangeFieldImp,m1*m2,n> {
  // Implement L : VD -> VR where VR is the space of FieldMatrix<RFI,m1,m2>
  // VD is space of FieldVector<DFT,n>
 public:
  //! \brief @copydoc FunctionSpaceBase::DomainFieldType 
  typedef DomainFieldImp DomainFieldType;
  //! \brief @copydoc FunctionSpaceBase::RangeFieldType 
  typedef RangeFieldImp RangeFieldType;
  //! \brief @copydoc FunctionSpaceBase::RangeType 
  typedef RangeMatrix<RangeFieldImp, m1,m2> RangeType;
  //! \brief type of base class 
  typedef FieldMatrix<RangeFieldImp,m1*m2,n> BaseType;
  //===== constructors
  /*! \brief Default constructor
   */
  MatrixMapping () : BaseType() {}
  /*! \brief Constructor initializing the whole matrix with a scalar
   */
  MatrixMapping (const RangeFieldImp& k) : BaseType(k) {}
  
  /** \brief returning reference to row 
      \param i number of row 
      \return Reference to row 
  */
  FieldVector<DomainFieldImp,n>& operator[](int i) {
    return static_cast<BaseType&>(*this)[i];
  }
  /** \brief returning reference to row 
      \param i number of row 
      \return Reference to row 
  */
  const FieldVector<DomainFieldImp,n>& operator[](int i) const {
    return static_cast<BaseType&>(*this)[i];
  }
};

//! \brief Traits class for matrix valued spaces 
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
struct MatrixSpaceTraits {
  //! \brief @copydoc FunctionSpaceBase::DomainFieldType 
  typedef DomainFieldImp DomainFieldType;
  //! \brief @copydoc FunctionSpaceBase::RangeFieldType 
  typedef RangeFieldImp RangeFieldType;
  //! \brief @copydoc FunctionSpaceBase::DomainType 
  typedef FieldVector<DomainFieldImp, n> DomainType;
  //! \brief @copydoc FunctionSpaceBase::RangeType 
  typedef RangeMatrix<RangeFieldImp, m1,m2> RangeType;
  //! \brief linear mapping type  
  typedef MatrixMapping<DomainFieldImp,RangeFieldImp, n, m1,m2> LinearMappingType;
  //! \brief scalar function space type 
  typedef MatrixFunctionSpace<DomainFieldImp,RangeFieldImp,n,1,1> ScalarFunctionSpaceType;
  //! \brief dimension of domain vector space 
  enum { DimRange = m1 * m2};
  //! \brief dimension of range vector space 
  enum { DimDomain = n };
};

//! \brief FunctionSpace class defining domain vector space and matrix valued range vector spaces */
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
class MatrixFunctionSpace : public
      FunctionSpaceBase<MatrixSpaceTraits<DomainFieldImp,RangeFieldImp,n,m1,m2> > {};


//! convert functions space to scalar function space
template <class FunctionSpaceImp>
struct ToScalarFunctionSpace {};

//! specialization for parameter list <domainfile,rangefield,dimDomain,dimRange> 
template <
  class DomainFieldImp, class RangeFieldImp, int dimDomain, int dimRange>
struct ToScalarFunctionSpace<
  FunctionSpace<DomainFieldImp, RangeFieldImp, dimDomain, dimRange> >
{
  typedef FunctionSpace<DomainFieldImp, RangeFieldImp, dimDomain, 1> Type;
};

} // end namespace Dune 
#endif
