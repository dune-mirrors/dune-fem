#ifndef DUNE_FUNCTIONSPACE_HH
#define DUNE_FUNCTIONSPACE_HH

#include <dune/fem/space/common/functionspaceinterface.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

// Forward declaration of
// base class for vector valued function spaces.
template< class DomainField, class RangeField, int dimD, int dimR >
class FunctionSpace;

//! \brief Traits class for vector function spaces 
template< class DomainField, class RangeField, int dimD,int dimR >
struct VectorSpaceTraits
{
  /** \copydoc Dune::FunctionSpaceInterface::DomainFieldType */
  typedef DomainField DomainFieldType;
  /** \copydoc Dune::FunctionSpaceInterface::RangeFieldType */
  typedef RangeField RangeFieldType;

  /** \brief dimension of range vector space */
  enum { dimDomain = dimD };
  /** \brief dimension of domain vector space */
  enum { dimRange = dimR };

  /** \copydoc Dune::FunctionSpaceInterface::DomainType */
  typedef FieldVector< DomainFieldType, dimDomain > DomainType;
  /** \copydoc Dune::FunctionSpaceInterface::RangeType */
  typedef FieldVector< RangeFieldType, dimRange> RangeType;
  
  /** \brief linear mapping type */
  typedef FieldMatrix< RangeFieldType, dimRange, dimDomain > LinearMappingType;

  /** \brief scalar function space type */
  typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >
    ScalarFunctionSpaceType;
};

/** @ingroup FunctionSpace
    \brief A vector valued function space.
   
    FunctionSpace defines what the types of the domain vector
    space and the range vector space for a function are. 
*/
template< class DomainField, class RangeField, int dimD, int dimR >
class FunctionSpace
: public FunctionSpaceInterface
    < VectorSpaceTraits< DomainField, RangeField, dimD, dimR > >
{
  typedef FunctionSpace< DomainField, RangeField, dimD, dimR > ThisType;

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
  typedef typename BaseType::row_type RowType;
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
    return static_cast<const BaseType&>(*this)[r][c];
  }
  /** \brief access to row r 
      \param[in] r row 
      \return reference to row r 
  */
  const RowType& row(int r) const{
    return static_cast<BaseType&>(*this)[r];
  }
 /** \brief access to row r 
      \param[in] r row 
      \return reference to row r 
 */
  RowType& row(int r){
    return static_cast<BaseType&>(*this)[r];
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
  /** \copydoc Dune::FunctionSpaceInterface::DomainFieldType */
  typedef DomainFieldImp DomainFieldType;
  /** \copydoc Dune::FunctionSpaceInterface::RangeFieldType */
  typedef RangeFieldImp RangeFieldType;
  /** \copydoc Dune::FunctionSpaceInterface::RangeType */
  typedef RangeMatrix<RangeFieldImp, m1,m2> RangeType;
  /** \brief type of base class */
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
    return static_cast<const BaseType&>(*this)[i];
  }
};

//! \brief Traits class for matrix valued spaces 
template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
struct MatrixSpaceTraits {
  /** \copydoc Dune::FunctionSpaceInterface::DomainFieldType */
  typedef DomainFieldImp DomainFieldType;
  /** \copydoc Dune::FunctionSpaceInterface::RangeFieldType */
  typedef RangeFieldImp RangeFieldType;
  /** \copydoc Dune::FunctionSpaceInterface::DomainType */
  typedef FieldVector<DomainFieldImp, n> DomainType;
  /** \copydoc Dune::FunctionSpaceInterface::RangeType */
  typedef RangeMatrix<RangeFieldImp, m1,m2> RangeType;
  /** \brief linear mapping type */
  typedef MatrixMapping<DomainFieldImp,RangeFieldImp, n, m1,m2> LinearMappingType;
  /** \brief scalar function space type */
  typedef MatrixFunctionSpace<DomainFieldImp,RangeFieldImp,n,1,1> ScalarFunctionSpaceType;
  /** \brief dimension of domain vector space */
  enum { dimRange = m1 * m2};
  /** \brief dimension of range vector space */
  enum { dimDomain = n };
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
