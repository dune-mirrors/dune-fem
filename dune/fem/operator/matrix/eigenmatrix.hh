#ifndef DUNE_FEM_EIGENMATRIX_HH
#define DUNE_FEM_EIGENMATRIX_HH

//- system includes
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>

//- local includes
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

    template <class T>
    class EigenMatrix
    {
      enum { defaultCol = -1 };
      enum { firstCol = defaultCol + 1 };

    public:
      typedef Eigen::SparseMatrix<T,Eigen::RowMajor> MatrixStorageType;
      typedef T Ttype;  //! remember the value type
      typedef EigenMatrix<T> ThisType;
      //! type of the base matrix (for consistency with ISTLMatrixObject)
      typedef ThisType MatrixBaseType;

    protected:
      MatrixStorageType matrix_;

    public:
      EigenMatrix(const EigenMatrix<T> &S) = delete;

      //! makes Matrix of zero length
      explicit EigenMatrix()
        : matrix_()
      {}

      //! make matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row
      //! and intialize all values with 'val' 
      EigenMatrix(int rows, int cols, int nz) // , const T& val = 0)
        : matrix_(rows,cols)
      {
        reserve(rows,cols,nz);
      }

      //! free memory for values_ and col_
      ~EigenMatrix()
      {}

      //! reserve memory for given rows, and number of non zeros,
      void reserve(int rows, int cols, int nz) 
      {
        matrix_.resize(rows,cols);
        matrix_.reserve(Eigen::VectorXi::Constant(cols,nz));;
      }

      //! return number of rows
      int rows() const
      {
        return matrix_.rows();
      }

      //! return number of columns
      int cols() const
      {
        return matrix_.cols();
      }
      
      //! set entry to value
      //! note, that every entry is performed into the matrix!
      //! also setting of value 0 will result in an entry. So these
      //! calls should be ommited on a higher level
      void set(int row, int col, T val)
      {
        matrix_.coeffRef(row,col) = val;
      }

      //! add value to row,col entry
      void add(int row, int col, T val)
      {
        matrix_.coeffRef(row,col) += val;
      }

      //! A(f) = ret, same as mult
      template <class DiscFType, class DiscFuncType>
      void apply(const DiscFType &f, DiscFuncType &ret) const
      {
        ret.dofVector().array().coefficients() = 
          matrix_ * f.dofVector().array().coefficients();
      }

      //! return value of entry (row,col)
      T operator() ( const int row, const int col ) const
      {
        T& val = const_cast<MatrixStorageType&>(matrix_).coeffRef(row,col);
        return val;
      }
      T operator() ( const unsigned int row, const unsigned int col ) const
      {
        return (*this)( int( row ), int( col ) );
      }
      T operator() ( const long unsigned int row, const long unsigned int col ) const
      {
        return this->operator()((unsigned int)(row), (unsigned int)(col) );
      }

      //! set all entries in row to zero
      void clearRow (int row)
      {
        std::cout << "EigenMatrix::clearRow not yet implemented" << std::endl;
        abort();
      }

      //! set all matrix entries to zero, no other value makes sense for
      //! sparse matrix
      void clear()
      {
        matrix_.setZero();
      }

      //! return max number of non zeros
      //! used in EigenMatrixObject::reserve
      int numNonZeros() const
      {
        return matrix_.nonZeros();
      }

      //! return number of non zeros in row
      //! used in ColCompMatrix::setMatrix
      // int numNonZeros(int i) const

      //! return pair< value, column >, used by BlockMatrix
      //! needed in ColCompMatrix::setMatrix
      // std::pair < const T , int > realValue(int index) const

      MatrixStorageType& data() 
      {
        return matrix_;
      }
      const MatrixStorageType& data() const 
      {
        return matrix_;
      }
    };

    template< class DomainSpace, class RangeSpace >
    struct EigenMatrixObject 
       : public SparseRowMatrixObject< DomainSpace, RangeSpace, EigenMatrix< typename DomainSpace :: RangeFieldType > >
    {
      typedef EigenMatrix< typename DomainSpace :: RangeFieldType > MatrixType;
      typedef SparseRowMatrixObject< DomainSpace, RangeSpace, MatrixType > BaseType;
      inline EigenMatrixObject( const DomainSpace &domainSpace,
                                const RangeSpace &rangeSpace,
                                const std::string &paramfile = "" )
        : BaseType( domainSpace, rangeSpace, paramfile )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPMATRIX_HH
