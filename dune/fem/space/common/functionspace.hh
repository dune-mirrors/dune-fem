#ifndef DUNE_FEM_FUNCTIONSPACE_HH
#define DUNE_FEM_FUNCTIONSPACE_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/dimension.hh>

#include <dune/fem/space/common/functionspaceinterface.hh>

namespace Dune
{

  namespace Fem
  {

    // Forward declaration of
    // base class for vector valued function spaces.
    template< class DomainField, class RangeField, int dimD, int dimR>
    class FunctionSpace;

    //! \brief Traits class for vector function spaces
    template< class DomainField, class RangeField, int dimD, int dimR>
    struct VectorSpaceTraits
    {
      /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainFieldType */
      typedef DomainField DomainFieldType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeFieldType */
      typedef RangeField RangeFieldType;

      /** \brief dimension of range vector space */
      enum { dimDomain = dimD };
      /** \brief dimension of domain vector space */
      enum { dimRange = dimR };

      /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainType */
      typedef FieldVector< DomainFieldType, dimDomain > DomainType;

      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeType */
      typedef FieldVector< RangeFieldType, dimRange> RangeType;

      /** \brief linear mapping type */
      typedef FieldMatrix< RangeFieldType, dimRange, dimDomain > LinearMappingType;

      /** \brief scalar function space type */
      typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1>
        ScalarFunctionSpaceType;
    };

    /** @ingroup FunctionSpace
        \brief A vector valued function space.

        FunctionSpace defines what the types of the domain vector
        space and the range vector space for a function are.
    */
    template< class DomainField, class RangeField, int dimD, int dimR>
    class FunctionSpace
    : public FunctionSpaceInterface
        < VectorSpaceTraits< DomainField, RangeField, dimD, dimR> >
    {
      typedef FunctionSpace< DomainField, RangeField, dimD, dimR> ThisType;

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
      /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainFieldType */
      typedef DomainFieldImp DomainFieldType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeFieldType */
      typedef RangeFieldImp RangeFieldType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeType */
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
      /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainFieldType */
      typedef DomainFieldImp DomainFieldType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeFieldType */
      typedef RangeFieldImp RangeFieldType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainType */
      typedef FieldVector<DomainFieldImp, n> DomainType;
      /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeType */
      typedef RangeMatrix<RangeFieldImp, m1,m2> RangeType;
      /** \brief linear mapping type */
      typedef MatrixMapping<DomainFieldImp,RangeFieldImp, n, m1,m2> LinearMappingType;
      /** \brief scalar function space type */
      typedef MatrixFunctionSpace<DomainFieldImp,RangeFieldImp,n,1,1> ScalarFunctionSpaceType;
      /** \brief dimension of range vector space */
      enum { dimRange = m1 * m2};
      /** \brief dimension of domain vector space */
      enum { dimDomain = n };
    };

    /*! @ingroup FunctionSpace
       \brief A matrix valued function space.
    */
    template <typename DomainFieldImp,typename RangeFieldImp,int n,int m1,int m2>
    class MatrixFunctionSpace : public
          FunctionSpaceInterface<MatrixSpaceTraits<DomainFieldImp,RangeFieldImp,n,m1,m2> > {};

    // function space converter objects
    // --------------------------------

    //! convert functions space to space with new dim domain
    template < class FunctionSpaceImp, int newDimDomain >
    struct ToNewDimDomainFunctionSpace;

    //! convert functions space to space with new dim range
    template < class FunctionSpaceImp, int newDimRange >
    struct ToNewDimRangeFunctionSpace;

    //! specialization for parameter list <domainfile,rangefield,dimDomain,dimRange,newDimDomain>
    template< class DomainFieldImp, class RangeFieldImp, int dimDomain, int dimRange, int newDimDomain >
    struct ToNewDimDomainFunctionSpace< FunctionSpace< DomainFieldImp, RangeFieldImp, dimDomain, dimRange >, newDimDomain >
    {
      typedef FunctionSpace< DomainFieldImp, RangeFieldImp, newDimDomain, dimRange> Type;
    };

    //! specialization for parameter list <domainfile,rangefield,dimDomain,dimRange,dimLocal>
    template< class DomainFieldImp, class RangeFieldImp, int n, int m1, int m2, int newDimDomain >
    struct ToNewDimDomainFunctionSpace< MatrixFunctionSpace< DomainFieldImp, RangeFieldImp, n, m1, m2 >, newDimDomain >
    {
      typedef MatrixFunctionSpace< DomainFieldImp, RangeFieldImp, newDimDomain, m1, m2 > Type;
    };

    //! specialization for parameter list <domainfile,rangefield,dimDomain,dimRange,dimLocal>
    template< class DomainFieldImp, class RangeFieldImp, int dimDomain, int dimRange, int newDimRange >
    struct ToNewDimRangeFunctionSpace< FunctionSpace< DomainFieldImp, RangeFieldImp, dimDomain, dimRange >, newDimRange >
    {
      typedef FunctionSpace< DomainFieldImp, RangeFieldImp, dimDomain, newDimRange > Type;
    };



    // GridFunctionSpace
    // -----------------

    namespace Impl
    {

      template< class GridPart, class >
      struct GridFunctionSpace;

      template< class GridPart, class K, int dimRange >
      struct GridFunctionSpace< GridPart, FunctionSpace< typename GridPart::ctype, K, GridPart::dimensionworld, dimRange > >
      {
        typedef FunctionSpace< typename GridPart::ctype, K, GridPart::dimensionworld, dimRange > Type;
      };

      template< class GridPart, class K, int rows, int cols >
      struct GridFunctionSpace< GridPart, MatrixFunctionSpace< typename GridPart::ctype, K, GridPart::dimensionworld, rows, cols > >
      {
        typedef MatrixFunctionSpace< typename GridPart::ctype, K, GridPart::dimensionworld, rows, cols > Type;
      };

      template< class GridPart, class K, int dimRange >
      struct GridFunctionSpace< GridPart, Dune::FieldVector< K, dimRange > >
      {
        typedef FunctionSpace< typename GridPart::ctype, K, GridPart::dimensionworld, dimRange > Type;
      };

      template< class GridPart, class K, int rows, int cols >
      struct GridFunctionSpace< GridPart, Dune::FieldMatrix< K, rows, cols > >
      {
        typedef MatrixFunctionSpace< typename GridPart::ctype, K, GridPart::dimensionworld, rows, cols > Type;
      };

      template< class GridPart, int dimRange >
      struct GridFunctionSpace< GridPart, Dune::Dim< dimRange > >
      {
        typedef typename GridFunctionSpace< GridPart, Dune::FieldVector< typename GridPart::ctype, dimRange > >::Type Type;
      };

    } // namespace Impl

    template< class GridPart, class T >
    using GridFunctionSpace = typename Impl::GridFunctionSpace< GridPart, T >::Type;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTIONSPACE_HH
