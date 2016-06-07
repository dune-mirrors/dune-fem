// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_COLCOMPSPMATRIX_HH
#define DUNE_COLCOMPSPMATRIX_HH

#if HAVE_DUNE_ISTL
#include <dune/istl/colcompmatrix.hh>
#else
namespace Dune
{
  template<class M>
  struct ColCompMatrix {};
} // namespace Dune
#endif

#include <dune/fem/operator/matrix/spmatrix.hh>

#include <vector>

namespace Dune
{
  /**
   *  @brief Converter for SparseRowMatrix to column-compressed matrix.
   *  Specialization for SparseRowMatrix
   */
  template<class B>
  class ColCompMatrix<Fem::SparseRowMatrix<B> >
  {
    public:
    /** @brief The type of the matrix converted. */
    typedef ColCompMatrix<Fem::SparseRowMatrix<B> > ThisType;
    /** @brief The type of the matrix to convert. */
    typedef Fem::SparseRowMatrix<B> Matrix;

    typedef typename Matrix::size_type size_type;

    /**
     * @brief Constructor that initializes the data.
     * @param mat The matrix to convert.
     */
    explicit ColCompMatrix(const Matrix& mat)
    {
      setMatrix(mat);
    }

    /** @brief Empty constructor. */
    ColCompMatrix() :
      N_(0), M_(0), Nnz_(0), values_(0), rowindex_(0), colstart_(0)
    {}

    /** @brief Destructor. */
    virtual ~ColCompMatrix()
    {
      if(N_+M_+Nnz_ != 0)
        free();
    }

    /** @brief Get the number of rows. */
    size_type N() const
    {
      return N_;
    }

    /** @brief Get the number of non zero entries. */
    size_type nnz() const
    {
      return Nnz_;
    }

    /** @brief Get the number of columns. */
    size_type M() const
    {
      return M_;
    }

    /** @brief Get the non-zero entries of the matrix. */
    B* getValues() const
    {
      return values_;
    }

    /** @brief Get the row indices of the non-zero entries of the matrix. */
    int* getRowIndex() const
    {
      return rowindex_;
    }

    /** @brief Get the column start indices. */
    int* getColStart() const
    {
      return colstart_;
    }

    ThisType& operator=(const Matrix& mat)
    {
      if(N_+M_+Nnz_ != 0)
        free();
      setMatrix(mat);
      return *this;
    }

    ThisType& operator=(const ThisType& mat)
    {
      if(N_+M_+Nnz_ != 0)
        free();
      N_=mat.N_;
      M_=mat.M_;
      Nnz_=mat.Nnz_;
      if(M_>0)
      {
        colstart_=new int[M_+1];
        for(size_type i=0; i<=M_; ++i)
          colstart_[i]=mat.colstart[i];
      }
      if(Nnz_>0)
      {
        values_ = new B[Nnz_];
        rowindex_ = new int[Nnz_];
        for(size_type i=0; i<Nnz_; ++i)
          values_[i]=mat.values[i];
        for(size_type i=0; i<Nnz_; ++i)
          rowindex_[i]=mat.rowindex[i];
      }
      return *this;
    }

    /** @brief Free allocated space. */
    void free()
    {
      delete[] values_;
      delete[] rowindex_;
      delete[] colstart_;
      N_=0;
      M_=0;
      Nnz_=0;
    }

    /** @brief Initialize data from given matrix. */
    virtual void setMatrix(const Matrix& mat)
    {
      N_=mat.rows();
      M_=mat.cols();

      // count the number of nonzeros per column
      colstart_=new int[M_+1];
      for(size_type i=0;i!=(M_+1);++i)
        colstart_[i]=0;
      Nnz_=0;
      size_type count(0);
      for(size_type i=0;i!=N_;++i)
      {
        Nnz_+=mat.numNonZeros(i);
        for(;count<(mat.numNonZeros()*(i+1));++count)
        {
          const auto pairIdx(mat.realValue(count));
          if(pairIdx.second!=Matrix::defaultCol)
            ++(colstart_[pairIdx.second+1]);
        }
      }

      // compute the starting positions and compute colstart
      std::vector<int> tempPos(M_,0);
      for(size_type i=1;i!=(M_+1);++i)
      {
        colstart_[i]+=colstart_[i-1];
        tempPos[i-1]=colstart_[i-1];
      }

      // fill the values and the index arrays
      values_=new B[Nnz_];
      rowindex_=new int[Nnz_];
      count=0;
      for(size_type i=0;i!=N_;++i)
      {
        for(size_type localCount=0;localCount<mat.numNonZeros();++localCount,++count)
        {
          const auto pairIdx(mat.realValue(count));
          if(pairIdx.second!=Matrix::defaultCol)
          {
            values_[tempPos[pairIdx.second]]=pairIdx.first;
            rowindex_[tempPos[pairIdx.second]]=i;
            ++(tempPos[pairIdx.second]);
          }
        }
      }

    }

    private:
    size_type N_;
    size_type M_;
    size_type Nnz_;
    B* values_;
    int* rowindex_;
    int* colstart_;
  };

}
#endif // #ifndef DUNE_COLCOMPSPMATRIX_HH

