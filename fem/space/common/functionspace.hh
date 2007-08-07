#ifndef DUNE_FUNCTIONSPACE_HH
#define DUNE_FUNCTIONSPACE_HH

#include <dune/fem/space/common/functionspaceinterface.hh>
#include <dune/common/fmatrix.hh>

namespace Dune{

// Forward declaration of
// base class for vector valued function spaces.
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
class FunctionSpace;

//! \brief Traits class for vector function spaces 
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m>
struct VectorSpaceTraits {
  //! \brief @copydoc FunctionSpaceInterface::DomainFieldType 
  typedef DomainFieldImp DomainFieldType;
  //! \brief @copydoc FunctionSpaceInterface::RangeFieldType 
  typedef RangeFieldImp RangeFieldType;
  //! \brief @copydoc FunctionSpaceInterface::DomainType 
  typedef FieldVector<DomainFieldImp, n> DomainType;
  //! \brief @copydoc FunctionSpaceInterface::RangeType 
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

/** @ingroup FunctionSpace
    \brief A vector valued function space.
   
    FunctionSpace defines what the types of the domain vector
    space and the range vector space for a function are. 
*/
template< typename DomainFieldImp, typename RangeFieldImp, int n, int m >
class FunctionSpace
: public FunctionSpaceInterface< VectorSpaceTraits< DomainFieldImp, RangeFieldImp, n, m > >
{
private:
  typedef FunctionSpace< DomainFieldImp, RangeFieldImp, n, m > ThisType;

public:
  typedef ThisType FunctionSpaceType;
};

/* Forward declaration of  
  base class for matrix valued function spaces.
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
  //! \brief @copydoc FunctionSpaceInterface::DomainFieldType 
  typedef DomainFieldImp DomainFieldType;
  //! \brief @copydoc FunctionSpaceInterface::RangeFieldType 
  typedef RangeFieldImp RangeFieldType;
  //! \brief @copydoc FunctionSpaceInterface::RangeType 
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
  //! \brief @copydoc FunctionSpaceInterface::DomainFieldType 
  typedef DomainFieldImp DomainFieldType;
  //! \brief @copydoc FunctionSpaceInterface::RangeFieldType 
  typedef RangeFieldImp RangeFieldType;
  //! \brief @copydoc FunctionSpaceInterface::DomainType 
  typedef FieldVector<DomainFieldImp, n> DomainType;
  //! \brief @copydoc FunctionSpaceInterface::RangeType 
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

/*! @ingroup FunctionSpace 
   \brief A matrix valued function space.
*/
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
class MatrixFunctionSpace : public
      FunctionSpaceInterface<MatrixSpaceTraits<DomainFieldImp,RangeFieldImp,n,m1,m2> > {};


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
