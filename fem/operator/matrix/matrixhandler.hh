#ifndef DUNE_MATRIXHANDLER_HH
#define DUNE_MATRIXHANDLER_HH

//- system includes 
#include <vector>

//- Dune istl includes 
#if HAVE_DUNE_ISTL
#include "istlmatrix.hh"
#endif

//- Dune includes 
#include <dune/fem/function/adaptivefunction.hh>

//- local includes 
#include "blockmatrix.hh"

namespace Dune {


  template <class SpaceType>
  class PreconditioningAdapter
  {
  public:  
    typedef AdaptiveDiscreteFunction<SpaceType> DiscreteFunctionType;
    typedef typename SpaceType :: RangeFieldType DofType;
  public:
    PreconditioningAdapter(std::string name, const SpaceType& space)
      : space_(space), diag_(name,space_) {} 
    
    //! apply this*x = ret 
    void precondition(const DofType * x, DofType * ret) const
    {
      const DofType* dofVec = diag_.leakPointer();
      const int vecsize = diag_.size();
      for(int i=0; i<vecsize; ++i) 
      {
        ret[i] = x[i]*dofVec[i];
      }
    } 
      
    //! multOEM interface 
    void multOEM (const DofType * x, DofType * ret) const  
    {
      precondition(x,ret);
    }   
      
    //! return false as this is left precondition
    bool rightPrecondition () const { return false; }

    //! clear matrix 
    void clear() { diag_.clear(); }

    DiscreteFunctionType& diag() { return diag_; }
  private:  
    const SpaceType& space_;
    DiscreteFunctionType diag_;

    PreconditioningAdapter(const PreconditioningAdapter&) {}
  };

  template <class SpaceType, class GradientSpaceType> 
  class MatrixHandlerSPMat
  {
    typedef typename SpaceType::GridType::template Codim<0>::Entity EntityType;
  public:  
    typedef SparseRowMatrix<double> MatrixType;
    typedef PreconditioningAdapter<SpaceType> PreconditionMatrixType;
    
    //! LocalMatrix 
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
      //! constructor taking entity and spaces for using mapToGlobal
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
        row_.resize(rowSpace.baseFunctionSet(rowEntity).numBaseFunctions());
        col_.resize(colSpace.baseFunctionSet(colEntity).numBaseFunctions());

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
      //! copy not allowed 
      MatrixNonSymetricHandle(const MatrixNonSymetricHandle &);

    public:
      //! return number of rows 
      int rows () const { return row_.size(); }
      //! return number of cols 
      int cols () const { return col_.size(); }

      //! add value to matrix entry
      void add(int localRow, int localCol , const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < (int) row_.size() );
        assert( localCol < (int) col_.size() );
        matrix_.add(row_[localRow],col_[localCol],value);
      }

      //! get matrix entry 
      double get(int localRow, int localCol) const 
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < (int) row_.size() );
        assert( localCol < (int) col_.size() );
        return matrix_(row_[localRow],col_[localCol]);
      }
      
      //! set matrix enrty to value 
      void set(int localRow, int localCol, const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < (int) row_.size() );
        assert( localCol < (int) col_.size() );
        matrix_.set(row_[localRow],col_[localCol],value);
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        const int row = rows();
        for(int i=0; i<row; ++i)
        {
          matrix_.clearRow( row_[i] );
        }
      }

      //! resort all global rows of matrix to have ascending numbering 
      void resort ()
      {
        const int row = rows();
        for(int i=0; i<row; ++i)
        {
          matrix_.resortRow( row_[i] );
        }
      }
    };

  public:
    typedef MatrixNonSymetricHandle<MatrixType> MatrixAddHandleType;

    const SpaceType & singleSpace_; 
    const GradientSpaceType & gradientSpace_;
    
    int singleMaxNumbers_;
    int gradMaxNumbers_;

    MatrixType stabMatrix_; 
    MatrixType gradMatrix_;
    MatrixType divMatrix_;
    MatrixType massMatrix_;
    PreconditionMatrixType* pcMatrix_;

    bool hasMassMatrix_;
    bool hasPcMatrix_;

    //! setup matrix handler 
    MatrixHandlerSPMat(const SpaceType & singleSpace, 
                       const GradientSpaceType & gradientSpace, 
                       const std::string& paramFile,
                       bool  hasMassMatrix = false)
      : singleSpace_(singleSpace)
      , gradientSpace_(gradientSpace) 
      , singleMaxNumbers_(-1)
      , gradMaxNumbers_ (-1) 
      , stabMatrix_()
      , gradMatrix_()
      , divMatrix_()
      , massMatrix_() // diagonal matrix 
      , pcMatrix_(0) 
      , hasMassMatrix_(hasMassMatrix)
      , hasPcMatrix_(false)
    {
      if(paramFile != "")
      {
        readParameter(paramFile,"Preconditioning",hasPcMatrix_);
      }
    }

    //! return reference to stability matrix 
    MatrixType & stabMatrix() { return stabMatrix_; }
    //! return reference to gradient matrix 
    MatrixType & gradMatrix() { return gradMatrix_; }
    //! return reference to divergence matrix 
    MatrixType & divMatrix()  { return divMatrix_; }
    //! return reference to mass matrix 
    MatrixType & massMatrix() { return massMatrix_; }
    //! return reference to preconditioning matrix 
    PreconditionMatrixType & pcMatrix() { 
      assert( pcMatrix_ );
      return *pcMatrix_; 
    }
    //! return true if mass matrix is used 
    bool hasMassMatrix() const { return hasMassMatrix_; }
    //! return true if preconditioning matrix is used 
    bool hasPcMatrix() const { return hasPcMatrix_; }

    //! resize all matrices and clear them 
    void clear() 
    {
      gradMatrix_.clear(); 
      divMatrix_.clear();
      stabMatrix_.clear();
      massMatrix_.clear();
      pcMatrix_.clear();
    }

    //! resize all matrices and clear them 
    void resize(bool verbose = false) 
    {
      if( ! hasBeenSetup() ) 
      {
        reserve(); 
      }
      else 
      {
        int singleSize = singleSpace_.size();
        int gradSize   = gradientSpace_.size();

        if(verbose)
          std::cout << "Resize Matrix with (" << singleSize << "," << gradSize << ")\n";

        gradMatrix_.resize(gradSize,singleSize);
        divMatrix_.resize(singleSize,gradSize);
        stabMatrix_.resize(singleSize,singleSize);

        if(hasMassMatrix())
        {
          massMatrix_.resize(gradSize,gradSize);
        }

        // pcMatrix_ is resized by dof manager 
        assert( (pcMatrix_) ? (pcMatrix_->size() == singleSize) : 1);
      }
    }

    //! clear mass matrix 
    void clearMass() 
    {
      if(hasMassMatrix())
      {
        massMatrix_.clear();
      }
    } 
    
    //! clear preconditioning matrix 
    void clearPcMatrix() 
    {
      if(hasPcMatrix())
      {
        assert( pcMatrix_ );
        pcMatrix_->clear();
      }
    } 

    //! returns true if memory has been reserved
    bool hasBeenSetup () const 
    {
      return (singleMaxNumbers_ > 0) && (gradMaxNumbers_ > 0);
    }

    //! reserve memory corresponnding to size of spaces 
    void reserve(bool verbose = false ) 
    {
      // if empty grid do nothing (can appear in parallel runs)
      if( (singleSpace_.begin()   != singleSpace_.end()) && 
          (gradientSpace_.begin() != gradientSpace_.end()) )
      {
        
        singleMaxNumbers_ = singleSpace_.baseFunctionSet(*(singleSpace_.begin())).numBaseFunctions();
        gradMaxNumbers_   = gradientSpace_.baseFunctionSet(*(gradientSpace_.begin())).numBaseFunctions();

        if(verbose) 
        {
          std::cout << "Reserve Matrix with (" << singleSpace_.size() << "," << gradientSpace_.size()<< ")\n";
          std::cout << "Number of base functions = (" << singleMaxNumbers_ << "," << gradMaxNumbers_ << ")\n";
        }

        assert( singleMaxNumbers_ > 0 );
        assert( gradMaxNumbers_ > 0 );

        // factor for non-conforming grid is 4 in 3d and 2 in 2d  
        //const int factor = (Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));
        const int factor = 1; //(Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));

        // upper estimate for number of neighbors 
        enum { dim = SpaceType :: GridType :: dimension };
        singleMaxNumbers_ *= (factor * dim * 2) + 1; // e.g. 7 for dim = 3
        gradMaxNumbers_   *= (factor * dim * 2) + 1;

        stabMatrix_.reserve(singleSpace_.size(),singleSpace_.size(),singleMaxNumbers_,0.0);
        gradMatrix_.reserve(gradientSpace_.size(),singleSpace_.size(),singleMaxNumbers_,0.0);
        divMatrix_.reserve(singleSpace_.size(),gradientSpace_.size(),gradMaxNumbers_,0.0);
        
        if( hasMassMatrix() ) 
        {
          massMatrix_.reserve(gradientSpace_.size(),gradientSpace_.size(),1,0.0);
        }

        if( hasPcMatrix() )
        {
          if(!pcMatrix_) pcMatrix_ = new PreconditionMatrixType("PCMatrix",singleSpace_); 
          pcMatrix_->clear();
        }
      }
    }

    //! resort row numbering in matrix to have ascending numbering 
    void resort() 
    {
      gradMatrix_.resort();
      divMatrix_.resort();
      stabMatrix_.resort();
      if( hasMassMatrix() )
      {
        massMatrix_.resort(); 
      }
    }

    void createPreconditionMatrix()
    {
      if(hasPcMatrix())
      {
        typedef typename PreconditionMatrixType :: DiscreteFunctionType
          DiscreteFunctionType; 
        DiscreteFunctionType & diag = pcMatrix().diag();

        if(hasMassMatrix())
        {
          divMatrix().getDiag( massMatrix(), gradMatrix() , diag );
        }
        else
        {
          divMatrix().getDiag( gradMatrix() , diag );
        }

        stabMatrix().addDiag( diag );

        double * diagPtr = diag.leakPointer();
        const int singleSize = singleSpace_.size();
        for(register int i=0; i<singleSize; ++i)
        {
          double val = diagPtr[i];
          // when using parallel Version , we could have zero on diagonal
          // for ghost elements 
          assert( (singleSpace_.grid().comm().size() > 1) ? 1 : (std::abs( val ) > 0.0 ) );
          if( std::abs( val ) > 0.0 )
          {
            val = 1.0/val;
            diagPtr[i] = val;
          }
        }
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
        if( matrix_.size(0) > 0 )
        {
          matrix_.add(rowIndex_,colIndex_, localMatrix_ );
        }
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

      //! set matrix enrty to value 
      void set(int localRow, int localCol, const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );

        assert( localRow < localMatrix_.rows() );
        assert( localCol < localMatrix_.cols() );

        localMatrix_[localRow][localCol] = value; 
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        localMatrix_.clear();
      }

      //! resort all global rows of matrix to have ascending numbering 
      void resort ()
      {
        //assert( matrix_.size(0) > 0 );
        //matrix_.resortRow( rowIndex_ );
      }
    };

  public:
    typedef MatrixNonSymetricHandle<MatrixType> MatrixAddHandleType;

    const SpaceType & singleSpace_; 
    const GradientSpaceType & gradientSpace_;
    int size_;

    int numSingleBaseFct_;
    int numGradBaseFct_;

    const int maxNumberUnknowns_;

    int singleMaxNumbers_;
    int gradMaxNumbers_;

    MatrixType * stabMatrix_; 
    MatrixType * divMatrix_;

    MatrixType gradMatrix_;
    MatrixType massMatrix_;

    MatrixType systemMatrix_; 

    bool hasMassMatrix_;
    bool hasPcMatrix_;

    MatrixHandlerBM(const SpaceType & singleSpace, 
                    const GradientSpaceType & gradientSpace,
                    bool hasMassMatrix, bool hasPcMatrix )
      : singleSpace_(singleSpace)
      , gradientSpace_(gradientSpace) 
      , size_(singleSpace_.indexSet().size(0))
      , maxNumberUnknowns_(10* (singleSpace_.baseFunctionSet(*(singleSpace_.begin())).numBaseFunctions()))
      , stabMatrix_(0)
      , divMatrix_ (0) 
      , hasMassMatrix_(hasMassMatrix)
      , hasPcMatrix_(hasPcMatrix)
    {
      reserve();
    }

    void reserve(bool verbose = false)
    {
      // if empty grid do nothing (can appear in parallel runs)
      if( (singleSpace_.begin()   != singleSpace_.end()) && 
          (gradientSpace_.begin() != gradientSpace_.end()) )
      {
        // get number of elements 
        size_ = singleSpace_.indexSet().size(0); 

        numSingleBaseFct_ = singleSpace_.baseFunctionSet(*(singleSpace_.begin())).numBaseFunctions();
        numGradBaseFct_   = gradientSpace_.baseFunctionSet(*(gradientSpace_.begin())).numBaseFunctions();

        if(verbose) 
        {
          std::cout << "Reserve Matrix with (" << size_ << "," << size_ << ")\n";
          std::cout << "Number of base functions = (" << numSingleBaseFct_ << "," << numGradBaseFct_ << ")\n";
        }

        singleMaxNumbers_ = numSingleBaseFct_;
        gradMaxNumbers_   = numGradBaseFct_;

        assert( singleMaxNumbers_ > 0 );
        assert( gradMaxNumbers_ > 0 );

        // factor for non-conforming grid is 4 in 3d and 2 in 2d  
        //const int factor = (Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));
        //const int factor = 1; //(Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));
        //(Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));

        // upper estimate for number of neighbors 
        enum { dim = SpaceType :: GridType :: dimension };

        // number of neighbors + 1 
        const int factor = (dim * 2) + 1; 
        // double stencil for system matrix 
        const int doubleFactor = (dim * 4) + 1; 
        
        singleMaxNumbers_ *= factor; // e.g. 7 for dim = 3
        gradMaxNumbers_   *= factor;

        if( !stabMatrix_ ) stabMatrix_ = new MatrixType ();
        stabMatrix().reserve(size_,size_,numSingleBaseFct_,numSingleBaseFct_,factor);
        
        if( !divMatrix_ ) divMatrix_ = new MatrixType ();
        divMatrix().reserve(size_,size_,numSingleBaseFct_,numGradBaseFct_,factor); 
        
        gradMatrix_.reserve(size_,size_,numGradBaseFct_, numSingleBaseFct_,factor); 

        systemMatrix_.reserve(size_,size_,numSingleBaseFct_,numSingleBaseFct_,2*doubleFactor); 

        if( hasMassMatrix() ) 
        {
          massMatrix_.reserve(size_,size_, numGradBaseFct_, numGradBaseFct_, 1 ); 
        }
      }
      else 
      {
        size_ = 0;
      }
    }

    MatrixType & stabMatrix() 
    { 
      assert( stabMatrix_ );
      return *stabMatrix_; 
    }
    
    MatrixType & gradMatrix() { return gradMatrix_; }
    MatrixType & divMatrix() 
    { 
      assert( divMatrix_ );
      return *divMatrix_; 
    }
    MatrixType & pcMatrix() 
    { 
      return systemMatrix_; 
    }
    MatrixType & massMatrix() { return massMatrix_; }

    MatrixType & systemMatrix() { return systemMatrix_; }

    //! return true if mass matrix is used 
    bool hasMassMatrix() const { return hasMassMatrix_; }
    //! return true if preconditioning matrix is used 
    bool hasPcMatrix() const { return hasPcMatrix_; }

    //! resize all matrices and clear them 
    void clear() 
    {
      gradMatrix_.clear(); 
      divMatrix().clear();
      stabMatrix().clear();
      massMatrix_.clear();
      systemMatrix_.clear();
    }

    void generateSystemMatrix () 
    {
      if( hasMassMatrix() )
      {
        MatrixType tmp;
        tmp.multiply( massMatrix_, gradMatrix_ );
        systemMatrix_.multiply( divMatrix() , tmp );
      }
      else 
      {
        systemMatrix_.multiply( divMatrix() , gradMatrix_ );
      }

      systemMatrix_.add( stabMatrix() );
      systemMatrix_.resort();

      delete divMatrix_; divMatrix_ = 0;
      delete stabMatrix_; stabMatrix_ = 0;
    }

    //! clear mass matrix 
    void clearMass() 
    {
      if(hasMassMatrix())
      {
        massMatrix_.clear();
      }
    } 
    
    //! clear preconditioning matrix 
    void clearPcMatrix() 
    {
    } 

    //! returns true if memory has been reserved
    bool hasBeenSetup () const 
    {
      return (singleMaxNumbers_ > 0) && (gradMaxNumbers_ > 0);
    }

    //! resize all matrices and clear them 
    void resize(bool verbose = false) 
    {
      if( ! hasBeenSetup() ) 
      {
        reserve(); 
      }
      else 
      {
        size_ = singleSpace_.indexSet().size(0); 
        if(verbose) 
        {
          std::cout << "Reserve Matrix with (" << size_ << "," << size_ << ")\n";
          std::cout << "Number of base functions = (" << numSingleBaseFct_ << "," << numGradBaseFct_ << ")\n";
        }

        gradMatrix_.resize(size_);

        if( !divMatrix_ ) divMatrix_ = new MatrixType ();
        divMatrix().resize(size_);

        if( !stabMatrix_ ) stabMatrix_ = new MatrixType ();
        stabMatrix().resize(size_);
      }
    }

    //! resort row numbering in matrix to have ascending numbering 
    void resort() 
    {
      gradMatrix_.resort();
      divMatrix().resort();
      stabMatrix().resort();
      if( hasMassMatrix() )
      {
        massMatrix_.resort(); 
      }
    }

    void createPreconditionMatrix() {}
  };

#if HAVE_DUNE_ISTL 
  template <class SpaceType, class GradientSpaceType> 
  class MatrixHandlerISTL 
  {
    typedef typename SpaceType::GridType::template Codim<0>::Entity EntityType;
  public:  

    //typedef FieldMatrix<double,1,1> LittleBlockType; 
    typedef FieldMatrix<double,6,6> LittleBlockType; 
    typedef ImprovedBCRSMatrix< LittleBlockType , SpaceType > MatrixType;
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
        row_.resize(rowSpace.baseFunctionSet(rowEntity).numBaseFunctions());
        col_.resize(colSpace.baseFunctionSet(colEntity).numBaseFunctions());

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
      , singleMax_(5 * (singleSpace_.baseFunctionSet(*(singleSpace_.begin())).numBaseFunctions()))
      , gradMax_(5 * (gradientSpace_.baseFunctionSet(*(singleSpace_.begin())).numBaseFunctions()))
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
      enum { dim = EntityType :: dimension };
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
