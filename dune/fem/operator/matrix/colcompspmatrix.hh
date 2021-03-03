// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_COLCOMPSPMATRIX_HH
#define DUNE_COLCOMPSPMATRIX_HH

#if HAVE_DUNE_ISTL
#include <dune/istl/bccsmatrixinitializer.hh>

// if the original deprecated header has not been included
// it is ok to declare the ColCompMatrix here
#ifndef DUNE_ISTL_COLCOMPMATRIX_HH
namespace Dune
{
  template <class M, class RowIndex = int>
  struct ColCompMatrix : public ISTL::Impl::BCCSMatrix< typename M::field_type, RowIndex>{};
} // namespace Dune
#endif

#else // #if HAVE_DUNE_ISTL
namespace Dune
{
  template<class M,class RowIndex=int>
  struct ColCompMatrix {};
} // namespace Dune
#endif // #if HAVE_DUNE_ISTL

#include <dune/fem/operator/matrix/spmatrix.hh>

#include <vector>

namespace Dune
{
  /**
   *  @brief Converter for SparseRowMatrix to column-compressed matrix.
   *  Specialization for SparseRowMatrix
   */
  template <class T, class IndexT,class ValuesVector, class IndicesVector,class RowIndex>
  class ColCompMatrix< Fem::SparseRowMatrix<T,IndexT,ValuesVector,IndicesVector>, RowIndex >
  {
    public:
    /** @brief The type of the matrix converted. */
    typedef ColCompMatrix< Fem::SparseRowMatrix<T,IndexT,ValuesVector,IndicesVector>> ThisType;
    /** @brief The type of the matrix to convert. */
    typedef Fem::SparseRowMatrix<T,IndexT,ValuesVector,IndicesVector> Matrix;

    typedef typename Matrix::size_type size_type;

    typedef RowIndex RowIndexType;

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
    T* getValues() const
    {
      return values_;
    }

    /** @brief Get the row indices of the non-zero entries of the matrix. */
    RowIndexType* getRowIndex() const
    {
      return rowindex_;
    }

    /** @brief Get the column start indices. */
    RowIndexType* getColStart() const
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
        colstart_=new RowIndexType[M_+1];
        for(size_type i=0; i<=M_; ++i)
          colstart_[i]=mat.colstart[i];
      }
      if(Nnz_>0)
      {
        values_ = new T[Nnz_];
        rowindex_ = new RowIndexType[Nnz_];
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
      colstart_= new RowIndexType[M_+1];
      for(size_type i=0;i!=(M_+1);++i)
        colstart_[i]=0;

      Nnz_ = 0;
      for(size_type row=0; row < N_; ++row)
      {
        const size_type endRow = mat.endRow( row );
        for( size_type col = mat.startRow( row ); col < endRow; ++col )
        {
          const auto pairIdx(mat.realValue( col ));
          if( pairIdx.second!=Matrix::defaultCol )
          {
            ++(colstart_[pairIdx.second+1]);
            ++Nnz_;
          }
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
      values_=new T[Nnz_];
      rowindex_=new RowIndexType[Nnz_];
      for(size_type row = 0; row < N_ ; ++ row)
      {
        const size_type endRow = mat.endRow( row );
        for( size_type col = mat.startRow( row ); col < endRow; ++col )
        {
          const auto pairIdx(mat.realValue( col ));
          if(pairIdx.second!=Matrix::defaultCol)
          {
            values_[tempPos[pairIdx.second]] = pairIdx.first;
            rowindex_[tempPos[pairIdx.second]] = row ;
            ++(tempPos[pairIdx.second]);
          }
        }
      }

    }

    private:
    size_type N_;
    size_type M_;
    size_type Nnz_;
    T* values_;
    RowIndexType* rowindex_;
    RowIndexType* colstart_;
  };

}
#endif // #ifndef DUNE_COLCOMPSPMATRIX_HH

