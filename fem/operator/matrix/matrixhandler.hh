#ifndef DUNE_MATRIXHANDLER_HH
#define DUNE_MATRIXHANDLER_HH

#include "blockmatrix.hh"

#if HAVE_DUNE_ISTL 
#include <dune/istl/bcrsmatrix.hh>
#endif

namespace Dune {

  template <class SpaceType, class GradientSpaceType> 
  class MatrixHandlerSPMat
  {
    typedef typename SpaceType::GridType::template Codim<0>::Entity EntityType;
  public:  
    typedef SparseRowMatrix<double> MatrixType;
    typedef MatrixType PreconditionMatrixType;
    
    template <class MatrixImp> 
    class MatrixNonSymetricHandle
    {
      typedef MatrixImp MatrixType;
      MatrixType & matrix_; 
      
      const int rowIndex_;
      const int colIndex_;

      std::vector<int> row_;
      std::vector<int> col_;
      
    public:  
      template <class EntityImp, class RowSpaceType, class ColSpaceType> 
      MatrixNonSymetricHandle(MatrixType & m,
                      const EntityImp & rowEntity,
                      const RowSpaceType & rowSpace,
                      const EntityImp & colEntity,
                      const ColSpaceType & colSpace)
        : matrix_(m)
        , rowIndex_(rowSpace.indexSet().index(rowEntity))
        , colIndex_(colSpace.indexSet().index(colEntity))
      {
        row_.resize(rowSpace.getBaseFunctionSet(rowEntity).numBaseFunctions());
        col_.resize(colSpace.getBaseFunctionSet(colEntity).numBaseFunctions());

        {
          const size_t rows = row_.size();
          for(size_t i=0; i<rows; ++i) 
            row_[i] = rowSpace.mapToGlobal( rowEntity, i );
        }
        {
          const size_t cols = col_.size();
          for(size_t i=0; i<cols; ++i) 
            col_[i] = colSpace.mapToGlobal( colEntity, i );
        }
      }

    private: 
      MatrixNonSymetricHandle(const MatrixNonSymetricHandle &);

    public:
      int rows () const { return row_.size(); }
      int cols () const { return col_.size(); }

      void add(int localRow, int localCol , const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < (int) row_.size() );
        assert( localCol < (int) col_.size() );
        matrix_.add(row_[localRow],col_[localCol],value);
      }

      double get(int localRow, int localCol) const 
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < (int) row_.size() );
        assert( localCol < (int) col_.size() );
        return matrix_(row_[localRow],col_[localCol]);
      }
    };

  public:
    typedef MatrixNonSymetricHandle<MatrixType> MatrixAddHandleType;

    const SpaceType & singleSpace_; 
    const GradientSpaceType & gradientSpace_;
    const int maxNumberUnknowns_;

    MatrixType stabMatrix_; 
    MatrixType gradMatrix_;
    MatrixType divMatrix_;
    MatrixType massMatrix_;
    MatrixType pcMatrix_;
    bool hasMassMatrix_;
    bool hasPcMatrix_;

    MatrixHandlerSPMat(const SpaceType & singleSpace, 
                       const GradientSpaceType & gradientSpace, 
                       bool  hasMassMatrix = false , 
                       bool hasPcMatrix = false )
      : singleSpace_(singleSpace)
      , gradientSpace_(gradientSpace) 
      , maxNumberUnknowns_(10* (singleSpace_.getBaseFunctionSet(*(singleSpace_.begin())).numBaseFunctions()))
      , stabMatrix_(singleSpace_.size(),singleSpace_.size(),maxNumberUnknowns_,0.0)
      , gradMatrix_(gradientSpace_.size(),singleSpace_.size(),maxNumberUnknowns_,0.0)
      , divMatrix_(singleSpace_.size(),gradientSpace_.size(),maxNumberUnknowns_,0.0)
      , massMatrix_(gradientSpace_.size(),gradientSpace_.size(),1,0.0) // diagonal matrix 
      , pcMatrix_(singleSpace_.size(),singleSpace_.size(),1,0.0) // diagonal matrix 
      , hasMassMatrix_(hasMassMatrix)
      , hasPcMatrix_(hasPcMatrix)
    {}

    MatrixType & stabMatrix() { return stabMatrix_; }
    MatrixType & gradMatrix() { return gradMatrix_; }
    MatrixType & divMatrix()  { return divMatrix_; }
    MatrixType & massMatrix() { return massMatrix_; }
    MatrixType & pcMatrix() { return pcMatrix_; }
    bool hasMassMatrix() const { return hasMassMatrix_; }
    bool hasPcMatrix() const { return hasPcMatrix_; }

    void resizeAndClear() 
    {
      int singleSize = singleSpace_.size();
      int gradSize   = gradientSpace_.size();
      std::cout << "Resize Matrix with " << singleSize << "\n";

      gradMatrix_.resize(gradSize,singleSize);
      gradMatrix_.clear();
        
      divMatrix_.resize(singleSize,gradSize);
      divMatrix_.clear();
  
      stabMatrix_.resize(singleSize,singleSize);
      stabMatrix_.clear();

      if(hasMassMatrix())
      {
        massMatrix_.resize(gradSize,gradSize);
        massMatrix_.clear();
      }

      if(hasPcMatrix()) 
      {
        pcMatrix_.resize(singleSize);
        pcMatrix_.clear();
      }
    }

    void clearMass() 
    {
      if(hasMassMatrix())
      {
        massMatrix_.clear();
      }
    } 
    
    void clearPcMatrix() 
    {
      if(hasPcMatrix())
      {
        pcMatrix_.clear();
      }
    } 
  };


  // --BlockMatrixHandle
  template <class SpaceType, class GradientSpaceType> 
  class MatrixHandlerBM 
  {
    typedef typename SpaceType::GridType::template Codim<0>::Entity EntityType;
  public:  
    typedef BlockMatrix<double> MatrixType;
    typedef MatrixType PreconditionMatrixType;
    
    template <class MatrixImp> 
    class MatrixNonSymetricHandle
    {
      typedef MatrixImp MatrixType;
      typedef DenseMatrix<double> LittleBlockType;

      MatrixType & matrix_; 
      
      const int rowIndex_;
      const int colIndex_;

      LittleBlockType & localMatrix_;

    public:  
      template <class EntityImp, class RowSpaceType, class ColSpaceType> 
      MatrixNonSymetricHandle(MatrixType & m,
                      const EntityImp & rowEntity,
                      const RowSpaceType & rowSpace,
                      const EntityImp & colEntity,
                      const ColSpaceType & colSpace)
        : matrix_(m)
        , rowIndex_(rowSpace.indexSet().index(rowEntity))
        , colIndex_(colSpace.indexSet().index(colEntity))
        , localMatrix_(matrix_.getMatrix())
      {
      }

      // add block to block matrix 
      ~MatrixNonSymetricHandle()
      {
        // finalize by adding block to global matrix 
        addLocalMatrix();
      }

    private: 
      MatrixNonSymetricHandle(const MatrixNonSymetricHandle &);

      void addLocalMatrix() 
      {
        matrix_.add(rowIndex_,colIndex_, localMatrix_ );
      }

    public:
      int rows () const { return localMatrix_.rows(); }
      int cols () const { return localMatrix_.cols(); }

      void add(int localRow, int localCol , const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < localMatrix_.rows() );
        assert( localCol < localMatrix_.cols() );

        localMatrix_[localRow][localCol] += value; 
      }

      double get(int localRow, int localCol ) const 
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < localMatrix_.rows() );
        assert( localCol < localMatrix_.cols() );

        return localMatrix_[localRow][localCol]; 
      }
    };

  public:
    typedef MatrixNonSymetricHandle<MatrixType> MatrixAddHandleType;

    const SpaceType & singleSpace_; 
    const GradientSpaceType & gradientSpace_;
    int size_;
    const int numSingleBaseFct_;
    const int numGradBaseFct_;

    const int maxNumberUnknowns_;

    MatrixType stabMatrix_; 
    MatrixType gradMatrix_;
    MatrixType divMatrix_;

    MatrixHandlerBM(const SpaceType & singleSpace, 
                    const GradientSpaceType & gradientSpace)
      : singleSpace_(singleSpace)
      , gradientSpace_(gradientSpace) 
      , size_(singleSpace_.indexSet().size(0))
      , numSingleBaseFct_(singleSpace_.getBaseFunctionSet( *singleSpace_.begin()).numBaseFunctions())
      , numGradBaseFct_(gradientSpace_.getBaseFunctionSet( *gradientSpace_.begin()).numBaseFunctions())
      , maxNumberUnknowns_(10* (singleSpace_.getBaseFunctionSet(*(singleSpace_.begin())).numBaseFunctions()))
      , stabMatrix_(size_,size_,
          numSingleBaseFct_,
          numSingleBaseFct_,
          maxNumberUnknowns_)
      , gradMatrix_(size_,size_
          , numGradBaseFct_
          , numSingleBaseFct_
          ,maxNumberUnknowns_)
      , divMatrix_(size_,size_
          , numSingleBaseFct_
          , numGradBaseFct_
          , maxNumberUnknowns_)
    {}

    MatrixType & stabMatrix() { return stabMatrix_; }
    MatrixType & gradMatrix() { return gradMatrix_; }
    MatrixType & divMatrix() { return divMatrix_; }

    void resizeAndClear() 
    {
      size_ = singleSpace_.indexSet().size(0); 

      gradMatrix_.resize(size_);
      gradMatrix_.clear();
        
      divMatrix_.resize(size_);
      divMatrix_.clear();
  
      stabMatrix_.resize(size_);
      stabMatrix_.clear();
    }
  };

#if HAVE_DUNE_ISTL 
  ///////////////////////////////////////////////////////
  // --BlockMatrixHandle
  //////////////////////////////////////////////////////
  template <class LittleBlockType> 
  class ImprovedBCRSMatrix : public BCRSMatrix<LittleBlockType> 
  {
      typedef BCRSMatrix<LittleBlockType> BaseType; 
      typedef typename BaseType :: RowIterator RowIteratorType ;
      typedef typename BaseType :: ColIterator ColIteratorType ;

      typedef typename BaseType :: size_type size_type;

      size_type nz_;

      int localRows_; 
      int localCols_;
    public:
      ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz)
        //: BaseType ()//rows,cols,rows*nz,bm)
        : BaseType (rows,cols, BaseType::row_wise)
        , nz_(nz)
      {
        std::cout << "Create Matrix with " << rows << " x " << cols << "\n";
      }

      void multOEM (const double * arg, double * dest)
      {
        RowIteratorType endi=this->end();
        for (RowIteratorType i=this->begin(); i!=endi; ++i)
        {
          int r = i.index();
          for(int k=0; k<localRows_; ++k)
          {
            int row = r * localRows_ + k; 
            dest[row] = 0.0;
          }

          ColIteratorType endj = (*i).end();
          for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
          {
            for(int k=0; k<localRows_; ++k) 
            {
              int row = r * localRows_ + k; 
              for(int c=0; c<localCols_; ++c)
              {
                int col = j.index() * localCols_ + c;  
                dest[row] += (*j)[k][c] * arg[col]; 
              }
            }
          }
        }
      }

      /*
      template <class RowSpaceType, class ColSpaceType> 
      void setup(const RowSpaceType & rowSpace, 
                 const ColSpaceType & colSpace) 
      {
        int rows = rowSpace.size();
        for(int i=0; i<rows; ++i) 
        {
          this->setrowsize(i,nz_);  
        }
        this->endrowsizes();

        typedef typename RowSpaceType :: IteratorType IteratorType; 
        typedef typename RowSpaceType :: GridType:: template Codim<0> ::
          Entity EntityType;
        
        IteratorType endit = rowSpace.end();
        for(IteratorType it = rowSpace.begin(); it != endit; ++it)
        {
          EntityType & en = *it;
          typedef typename RowSpaceType :: BaseFunctionSetType RowBaseSetType;
          typedef typename ColSpaceType :: BaseFunctionSetType ColBaseSetType;
          
          const RowBaseSetType & rowBase = rowSpace.getBaseFunctionSet(en);
          const ColBaseSetType & colBase = colSpace.getBaseFunctionSet(en);

          const int localRows = rowBase.numBaseFunctions();
          const int localCols = colBase.numBaseFunctions();

          for(int row=0 ; row<localRows; ++row) 
          {
            int globalRow = rowSpace.mapToGlobal(en,row);
            for(int col=0; col<localCols; ++col)
            {
              int globalCol = colSpace.mapToGlobal(en,col);
              this->addindex( globalRow , globalCol );  
            }
          }
        }
        this->endindices();
      }
      */
      template <class RowSpaceType, class ColSpaceType> 
      void setup(const RowSpaceType & rowSpace, 
                 const ColSpaceType & colSpace) 
      {
        {
        typedef typename BaseType :: CreateIterator CreateIteratorType; 

        CreateIteratorType create = this->createbegin();
        CreateIteratorType endcreate = this->createend();

        typedef typename RowSpaceType :: IteratorType IteratorType; 
        typedef typename RowSpaceType :: GridType:: template Codim<0> ::
          Entity EntityType;
        typedef typename RowSpaceType :: GridPartType:: IntersectionIteratorType 
          IntersectionIteratorType; 
        
        IteratorType endit = rowSpace.end();
        for(IteratorType it = rowSpace.begin(); it != endit; ++it)
        {
          assert( create != endcreate );

          EntityType & en = *it;

          localRows_ = rowSpace.getBaseFunctionSet(en).numBaseFunctions();
          localCols_ = colSpace.getBaseFunctionSet(en).numBaseFunctions();

          const int elIndex = rowSpace.indexSet().index(en);
          create.insert( elIndex );

          typedef typename RowSpaceType :: BaseFunctionSetType RowBaseSetType;
          typedef typename ColSpaceType :: BaseFunctionSetType ColBaseSetType;

          IntersectionIteratorType endnit = rowSpace.gridPart().iend(en);
          for(IntersectionIteratorType nit =
              rowSpace.gridPart().ibegin(en); nit != endnit; ++nit)
          {
            if(nit.neighbor())
            {
              const int nbIndex = rowSpace.indexSet().index( *nit.outside() );

              create.insert( nbIndex );
            }
          }
          ++create;
        }

      }

        {
          RowIteratorType endi=this->end();
          for (RowIteratorType i=this->begin(); i!=endi; ++i)
          {
            ColIteratorType endj = (*i).end();
            for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
            {
              (*j) = 0.0;
            }
          }
        }
      }
  };

  template <class SpaceType, class GradientSpaceType> 
  class MatrixHandlerISTL 
  {
    typedef typename SpaceType::GridType::template Codim<0>::Entity EntityType;
  public:  

    //typedef FieldMatrix<double,1,1> LittleBlockType; 
    typedef FieldMatrix<double,6,6> LittleBlockType; 
    typedef ImprovedBCRSMatrix< LittleBlockType > MatrixType;
    typedef MatrixType PreconditionMatrixType;
    
    template <class MatrixImp> 
    class MatrixNonSymetricHandle
    {
      typedef MatrixImp MatrixType;
      typedef typename MatrixType:: block_type LittleBlockType;
      
      //MatrixType & matrix_; 
      
      //std::vector<int> row_;
      //std::vector<int> col_;

      //int rows_;
      //int cols_;
      int rowIndex_;
      int colIndex_;

      LittleBlockType & matrix_;
      
    public:  
      template <class EntityImp, class RowSpaceType, class ColSpaceType> 
      MatrixNonSymetricHandle(MatrixType & m,
                      const EntityImp & rowEntity,
                      const RowSpaceType & rowSpace,
                      const EntityImp & colEntity,
                      const ColSpaceType & colSpace)
        //: matrix_(m)
        : rowIndex_(rowSpace.indexSet().index(rowEntity))
        , colIndex_(colSpace.indexSet().index(colEntity))
        , matrix_(m[rowIndex_][colIndex_])
      {
        /*
        row_.resize(rowSpace.getBaseFunctionSet(rowEntity).numBaseFunctions());
        col_.resize(colSpace.getBaseFunctionSet(colEntity).numBaseFunctions());

        {
          const size_t rows = row_.size();
          for(size_t i=0; i<rows; ++i) 
            row_[i] = rowSpace.mapToGlobal( rowEntity, i );
        }
        {
          const size_t cols = col_.size();
          for(size_t i=0; i<cols; ++i) 
            col_[i] = colSpace.mapToGlobal( colEntity, i );
        }
        */
      }

    private: 
      MatrixNonSymetricHandle(const MatrixNonSymetricHandle &);

    public:
      void add(int localRow, int localCol , const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        /*
        assert( localRow < (int) row_.size() );
        assert( localCol < (int) col_.size() );

        const int row = row_[localRow];
        const int col = col_[localCol];

        assert( row < rows_ );
        assert( col < cols_ );
        */

        //matrix_.addindex( row , col ) ;
        //matrix_[row][col] += value;
        matrix_[localRow][localCol] += value;
      }
    };

  public:
    typedef MatrixNonSymetricHandle<MatrixType> MatrixAddHandleType;

    const SpaceType & singleSpace_; 
    const GradientSpaceType & gradientSpace_;

    const int singleMax_;
    const int gradMax_;

    MatrixType stabMatrix_; 
    MatrixType gradMatrix_;
    MatrixType divMatrix_;

    MatrixHandlerISTL(const SpaceType & singleSpace, 
                      const GradientSpaceType & gradientSpace)
      : singleSpace_(singleSpace)
      , gradientSpace_(gradientSpace) 
      , singleMax_(5 * (singleSpace_.getBaseFunctionSet(*(singleSpace_.begin())).numBaseFunctions()))
      , gradMax_(5 * (gradientSpace_.getBaseFunctionSet(*(singleSpace_.begin())).numBaseFunctions()))
      //, stabMatrix_(singleSpace_.size(),singleSpace_.size(), singleMax_)
      //, gradMatrix_(gradientSpace_.size(),singleSpace_.size(),singleMax_)
      //, divMatrix_(singleSpace_.size(),gradientSpace_.size(), gradMax_)
      , stabMatrix_(singleSpace_.indexSet().size(0),singleSpace_.indexSet().size(0), singleMax_)
      , gradMatrix_(gradientSpace_.indexSet().size(0),singleSpace_.indexSet().size(0),singleMax_)
      , divMatrix_(singleSpace_.indexSet().size(0),gradientSpace_.indexSet().size(0), gradMax_)
    {
      stabMatrix_.setup(singleSpace_,singleSpace_); 
      gradMatrix_.setup(gradientSpace_,singleSpace_); 
      divMatrix_.setup(singleSpace_,gradientSpace_); 
    }

    MatrixType & stabMatrix() { return stabMatrix_; }
    MatrixType & gradMatrix() { return gradMatrix_; }
    MatrixType & divMatrix() { return divMatrix_; }

    void resizeAndClear() 
    {
      {
        MatrixType m(singleSpace_.indexSet().size(0),singleSpace_.indexSet().size(0), singleMax_);
        m.setup(singleSpace_,singleSpace_);
        stabMatrix_= m;
      }
      
      { 
        MatrixType g (gradientSpace_.indexSet().size(0),singleSpace_.indexSet().size(0),singleMax_);
        g.setup(gradientSpace_,singleSpace_);
        gradMatrix_ = g;
      }

      {
        MatrixType d (singleSpace_.indexSet().size(0),gradientSpace_.indexSet().size(0), gradMax_);
        d.setup(singleSpace_,gradientSpace_);
        divMatrix_ = d;
      }
    }
  };
#endif // HAVE_DUNE_ISTL

  template <class MatrixImp,class LocalMatrixType, class RowSpaceImp, class ColSpaceImp>
  class MatrixDataHandler 
   : public CommDataHandleIF< 
      MatrixDataHandler<MatrixImp,LocalMatrixType,RowSpaceImp,ColSpaceImp>,
      typename RowSpaceImp::RangeFieldType > 
  {
  public:  
    typedef MatrixImp MatrixType;
    typedef typename RowSpaceImp::RangeFieldType DataType;
    typedef RowSpaceImp RowSpaceType;
    typedef ColSpaceImp ColSpaceType;

    typedef typename RowSpaceType :: GridType GridType;

  private:  
    
    MatrixType & matrix_;
    const RowSpaceImp & rowSpace_;
    const ColSpaceImp & colSpace_; 
  public:
    MatrixDataHandler(MatrixType & matrix, 
                      const RowSpaceImp & rowSpace,
                      const ColSpaceImp & colSpace)
      : matrix_(matrix)
      , rowSpace_(rowSpace)
      , colSpace_(colSpace)
    {
    }
    
    MatrixDataHandler(const MatrixDataHandler & org) 
      : matrix_(org.matrix_)
      , rowSpace_(org.rowSpace_)
      , colSpace_(org.colSpace_)
    {
    }

    bool contains (int dim, int codim) const
    {
      return (codim == 0);
    }

    bool fixedsize (int dim, int codim) const
    {
      return true;
    }

    template<class MessageBufferImp, class EntityType>
    void gather (MessageBufferImp& buff, const EntityType& en) const
    {
      {
        LocalMatrixType localMat(matrix_,en, rowSpace_, en, colSpace_ );
        writeToBuff(buff,localMat);
      }

      typedef typename RowSpaceType :: GridPartType GridPartType; 
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType; 
     
      const GridPartType & gridPart = rowSpace_.gridPart();
      
      {
        IntersectionIteratorType endnit = gridPart.iend(en);
        for (IntersectionIteratorType nit = gridPart.ibegin(en); nit != endnit; ++nit)
        {
          if(nit.neighbor())
          {
            typedef typename GridPartType :: GridType :: template
              Codim<0>::EntityPointer EntityPointerType; 
            typedef typename GridPartType :: GridType :: template
              Codim<0>::Entity EntityImp; 
            EntityPointerType neighEp = nit.outside();
            EntityImp&            nb = *neighEp;

            LocalMatrixType localMat(matrix_,en, rowSpace_, nb, colSpace_ );
            writeToBuff(buff,localMat);
          }
        }
      }
    }

    template<class MessageBufferImp, class EntityType>
    void scatter (MessageBufferImp& buff, const EntityType& en, size_t n)
    {
      {
        LocalMatrixType localMat(matrix_,en, rowSpace_, en, colSpace_ );
        readFromBuff(buff,localMat);
      }
      
      typedef typename RowSpaceType :: GridPartType GridPartType; 
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType; 
     
      const GridPartType & gridPart = rowSpace_.gridPart();
     
      {
        IntersectionIteratorType endnit = gridPart.iend(en);
        for (IntersectionIteratorType nit = gridPart.ibegin(en); nit != endnit; ++nit)
        {
          if(nit.neighbor())
          {
            typedef typename GridPartType :: GridType :: template
              Codim<0>::EntityPointer EntityPointerType; 
            typedef typename GridPartType :: GridType :: template
              Codim<0>::Entity EntityImp; 
            EntityPointerType neighEp = nit.outside();
            EntityImp&  nb = *neighEp;
            
            LocalMatrixType localMat(matrix_,en, rowSpace_, nb, colSpace_ );
            readFromBuff(buff,localMat);
          }
        }
      }
    }

    template<class EntityType>
    size_t size (const EntityType& en) const
    {
      LocalMatrixType localMat(matrix_,en, rowSpace_, en, colSpace_ );
      enum { factor = ((dim)*2 + 1) };
      size_t s = localMat.rows() * localMat.cols();
      s *= factor;
      return s;
    }
  private:  
    template <class MessageBufferImp>
    void writeToBuff(MessageBufferImp& buff, const LocalMatrixType &
        localMat ) const
    {
      const int rows = localMat.rows();
      const int cols = localMat.cols();
      for(int i=0; i<rows; ++i) 
      {
        for(int j=0; j<cols; ++j)
        {
          buff.write( localMat.get(i,j) );
        }
      }
    }
    
    template <class MessageBufferImp>
    void readFromBuff(MessageBufferImp& buff, LocalMatrixType &
        localMat ) 
    {
      const int rows = localMat.rows();
      const int cols = localMat.cols();
      DataType val;
      for(int i=0; i<rows; ++i) 
      {
        for(int j=0; j<cols; ++j)
        {
          buff.read( val );
          localMat.add(i,j,val);
        }
      }
    }
  };

  class EmptyDataHandler 
  {
  public:
    EmptyDataHandler() {}
    EmptyDataHandler(const EmptyDataHandler &) {}

    static EmptyDataHandler & instance () 
    {
      static EmptyDataHandler emptyData;
      return emptyData;
    }

    bool contains (int dim, int codim) const
    {
      return false;
    }

    bool fixedsize (int dim, int codim) const
    {
      return true;
    }

    template<class MessageBufferImp, class EntityType>
    void gather (MessageBufferImp& buff, const EntityType& en) const
    {
    }

    template<class MessageBufferImp, class EntityType>
    void scatter (MessageBufferImp& buff, const EntityType& en, size_t n)
    {
    }

    template<class EntityType>
    size_t size (const EntityType& en) const
    {
      return 0;
    }
  };

} // end namespace Dune
#endif
