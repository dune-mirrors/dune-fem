#ifndef DUNE_FEM_BLOCKSPMATRIX_HH
#define DUNE_FEM_BLOCKSPMATRIX_HH

//- system includes
#include <stack>

//- local includes
#include "spmatrix.hh"

namespace Dune
{

  namespace Fem
  {

    //*****************************************************************
    //
    //  --DenseMatrix
    //
    //! \brief DenseMatrix based on std::vector< std::vector< T > >
    //*****************************************************************
    template <class T>
    class DenseMatrix
    {
    public:
      typedef T Ttype;  //! remember the value type

      typedef std::vector < T > RowType;
    private:
      typedef std::vector < RowType > MatrixType;

      MatrixType matrix_;

      int rows_;
      int cols_;

    public:
      // creating empty matrix
      DenseMatrix() : //! makes Matrix of zero length
        matrix_()
      {
        rows_ = 0;
        cols_ = 0;
      }

      //! Copy Constructor
      DenseMatrix(const DenseMatrix<T> & org) :
        matrix_(org.matrix_) , rows_(org.rows_) , cols_(org.cols_)
      {
      }

      //! make matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row
      //! and intialize all values with 'val'
      DenseMatrix(int rows, int cols)
      {
        resize(rows,cols);
      }

      void resize(int rows, int cols)
      {
        if( (rows == rows_) && (cols == cols_) ) return ;

        matrix_.resize(rows);
        rows_ = rows;
        cols_ = cols;

        for(int row=0; row<rows_; ++row)
        {
          matrix_[row].resize(cols);
          for(int col=0; col<cols_; ++col) matrix_[row][col] = 0;
        }
      }

    /*******************************/
    /*  Access and info functions  */
    /*******************************/
      int rows() const {return rows_; }
      int cols() const {return cols_; }

      //int base() const {return base_;}
      T& operator() (int row, int col)
      {
        assert( row>= 0);
        assert( col>= 0);

        assert( row < rows_);
        assert( col < cols_);
        return matrix_[row][col];
      }
      const T& operator() (int row, int col) const
      {
        assert( row >= 0);
        assert( col >= 0);

        assert( row < rows_);
        assert( col < cols_);
        return matrix_[row][col];
      }

      RowType & operator [] (int row) { return matrix_[row]; }
      const RowType & operator [] (int row) const { return matrix_[row]; }

      // result = this * vec
      void mult(const T * vec, RowType & result) const
      {
        assert( ((int) result.size() != rows()) ?
            (std::cout << result.size() << " s|r " << rows() << "\n",0) : 1 );
        const int nRow= rows();
        const int nCol= cols();
        for(int row=0; row<nRow; ++row)
        {
          result[row] = 0;
          for(int col=0; col<nCol; ++col)
          {
            result[row] += matrix_[row][col]*vec[col];
          }
        }
      }

      // result = this * vec
      void multOEM(const T * vec, T * result) const
      {
        const int nRow= rows();
        const int nCol= cols();
        for(int row=0; row<nRow; ++row)
        {
          result[row] = 0;
          for(int col=0; col<nCol; ++col)
          {
            result[row] += matrix_[row][col]*vec[col];
          }
        }
      }

      // result = this * vec
      void mult(const RowType & vec, RowType & result) const
      {
        this->mult(&vec[0],result);
      }

      // result = this^T * vec
      void multTransposed(const RowType & vec, RowType & result) const
      {
        assert( (int) result.size() == cols() );
        const int nCols = cols();
        const int nRows = rows();
        for(int col=0; col<nCols; ++col)
        {
          result[col] = 0;
          for(int row=0; row<nRows; ++row)
          {
            result[col] += matrix_[row][col]*vec[row];
          }
        }
      }

      // this = A * B
      void multiply(const DenseMatrix & A, const DenseMatrix & B)
      {
        assert( A.cols() == B.rows() );

        resize( A.rows() , B.cols() );

        const int nRows = rows();
        const int nCols = cols();
        const int Acols = A.cols();
        for(int row=0; row<nRows; ++row)
        {
          for(int col=0; col<nCols; ++col)
          {
            T sum = 0;
            for(int k=0; k<Acols; ++k)
            {
              sum += A[row][k] * B[k][col];
            }
            matrix_[row][col] = sum;
          }
        }
      };

      // this = A * B
      void multiplyTransposed(const DenseMatrix & A, const DenseMatrix & B)
      {
        assert( A.cols() == B.cols() );

        resize(A.rows() , B.rows());

        for(int row=0; row<rows(); ++row)
        {
          for(int col=0; col<cols(); ++col)
          {
            T sum = 0;
            for(int k=0; k<A.cols(); ++k)
            {
              sum += A[row][k] * B[col][k];
            }
            matrix_[row][col] = sum;
          }
        }
      };

      //! this = A^T * A
      void multiply_AT_A(const DenseMatrix & A)
      {
        resize(A.cols() , A.cols());

        for(int row=0; row<rows(); ++row)
        {
          for(int col=0; col<cols(); ++col)
          {
            T sum = 0;
            const int aRows = A.rows();
            for(int k=0; k<aRows; ++k)
            {
              sum += A[k][row] * A[k][col];
            }
            matrix_[row][col] = sum;
          }
        }
      }

      //! scale matrix with scalar
      void scale ( const T& val )
      {
        for(int row=0; row<rows(); ++row)
        {
          for(int col=0; col<cols(); ++col)
          {
            matrix_[row][col] *= val;
          }
        }
      }

      //! set all values of the matrix to given value
      DenseMatrix<T> & operator = (const T & val)
      {
        for(int row=0; row<rows(); ++row)
        {
          for(int col=0; col<cols(); ++col)
          {
            matrix_[row][col] = val;
          }
        }
        return *this;
      }

      //! set all values of the matrix to values of given matrix
      DenseMatrix<T> & operator = (const DenseMatrix & org)
      {
        rows_ = org.rows_;
        cols_ = org.cols_;

        matrix_ = org.matrix_;
        return *this;
      }

      //! add matrix
      DenseMatrix<T> & operator += (const DenseMatrix & org)
      {
        const int nRows = rows();
        const int nCols = cols();
        assert( nRows == org.rows() );
        assert( nCols == org.cols() );
        for(int row=0; row<nRows; ++row)
        {
          for(int col=0; col<nCols; ++col)
          {
            matrix_[row][col] += org.matrix_[row][col];
          }
        }
        return *this;
      }

      //! substract matrix
      DenseMatrix<T> & operator -= (const DenseMatrix & org)
      {
        assert( rows() == rows() );
        assert( cols() == org.cols() );
        for(int row=0; row<rows(); ++row)
        {
          for(int col=0; col<cols(); ++col)
            matrix_[row][col] -= org.matrix_[row][col];
        }
        return *this;
      }

      //! print matrix
      void print(std::ostream & s=std::cout) const
      {
        for(int i=0; i<rows(); ++i)
        {
          for(int j=0; j<cols(); ++j)
            s << matrix_[i][j] << " ";
          s << std::endl;
        }
      }

      // set all matrix entries to zero
      void clear()
      {
        for(int row=0; row<rows(); ++row)
        {
          for(int col=0; col<cols(); ++col)
            matrix_[row][col] = 0;
        }
      }
    }; // end class DenseMatrix

    //! Send vector to output stream
    template<typename K>
    std::ostream& operator<< (std::ostream& s, const DenseMatrix<K> & matrix)
    {
      matrix.print(s);
      return s;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BLOCKSPMATRIX_HH
