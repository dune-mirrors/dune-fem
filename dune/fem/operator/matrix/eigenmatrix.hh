#ifndef DUNE_FEM_EIGENMATRIX_HH
#define DUNE_FEM_EIGENMATRIX_HH

#ifdef HAVE_EIGEN

// system includes
#include <iostream>
#include <string>
#include <utility>

// local includes
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/storage/objectstack.hh>
#include <Eigen/Sparse>

namespace Dune
{

  namespace Fem
  {

    //! EigenMatrix
    template <class T>
    class EigenMatrix
    {
      static constexpr int defaultCol = -1;
      static constexpr int firstCol = defaultCol + 1;

    public:
      //! matrix field type
      typedef T field_type;
      //! matrix index type
      typedef int size_type;
      typedef Eigen::SparseMatrix<field_type,Eigen::RowMajor> MatrixStorageType;
      typedef EigenMatrix<field_type,size_type> ThisType;
      //! type of the base matrix
      //! for consistency with ISTLMatrixObject
      typedef ThisType MatrixBaseType;

      EigenMatrix(const ThisType& ) = delete;

      //! construct matrix of zero size
      explicit EigenMatrix() :
        matrix_()
      {}

      //! construct matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row
      EigenMatrix(size_type rows, size_type cols, size_type nz) :
        matrix_(rows,cols)
      {
        reserve(rows,cols,nz);
      }

      //! reserve memory for given rows, columns and number of non zeros
      void reserve(size_type rows, size_type cols, size_type nz)
      {
        matrix_.resize(rows,cols);
        matrix_.reserve(Eigen::VectorXi::Constant(cols,nz));;
      }

      //! return number of rows
      size_type rows() const
      {
        return matrix_.rows();
      }

      //! return number of columns
      size_type cols() const
      {
        return matrix_.cols();
      }

      //! set entry to value (also setting 0 will result in an entry)
      void set(size_type row, size_type col, field_type val)
      {
        matrix_.coeffRef(row,col) = val;
      }

      //! add value to row,col entry
      void add(size_type row, size_type col, field_type val)
      {
        matrix_.coeffRef(row,col) += val;
      }

      //! ret = A*f
      template<class ArgDFType, class DestDFType>
      void apply(const ArgDFType& f, DestDFType& ret) const
      {
        ret.dofVector().array().coefficients() =
          matrix_ * f.dofVector().array().coefficients();
      }

      //! return value of entry (row,col)
      field_type operator()(size_type row, size_type col) const
      {
        return matrix_.coeffRef(row,col);
      }

      //! set all matrix entries to zero
      void clear()
      {
        matrix_.setZero();
      }

      //! set all entries in row to zero
      void clearRow (size_type row)
      {
        std::cout << "EigenMatrix::clearRow not yet implemented" << std::endl;
        abort();
      }

      //! return max number of non zeros
      //! used in EigenMatrixObject::reserve
      size_type numNonZeros() const
      {
        return matrix_.nonZeros();
      }

      //! return number of non zeros in row
      //! used in ColCompMatrix::setMatrix
      size_type numNonZeros(size_type i) const
      {
        std::cout << "EigenMatrix::numNonZeros not yet implemented" << std::endl;
        abort();
      }

      //! return pair (value,column)
      //! used in ColCompMatrix::setMatrix
      std::pair<const field_type, size_type> realValue(size_type index) const
      {
        std::cout << "EigenMatrix::realValue not yet implemented" << std::endl;
        abort();
      }

      MatrixStorageType& data()
      {
        return matrix_;
      }

      const MatrixStorageType& data() const
      {
        return matrix_;
      }

    protected:
      MatrixStorageType matrix_;
    };

    template< class DomainSpace, class RangeSpace >
    struct EigenMatrixObject
       : public SparseRowMatrixObject< DomainSpace, RangeSpace, EigenMatrix< typename DomainSpace::RangeFieldType> >
    {
      typedef EigenMatrix< typename DomainSpace::RangeFieldType > MatrixType;
      typedef SparseRowMatrixObject< DomainSpace, RangeSpace, MatrixType > BaseType;
      inline EigenMatrixObject( const DomainSpace &domainSpace,
                                const RangeSpace &rangeSpace,
                                const std::string &paramfile = "" )
        : BaseType( domainSpace, rangeSpace, paramfile )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif

#endif // #ifndef DUNE_FEM_SPMATRIX_HH
