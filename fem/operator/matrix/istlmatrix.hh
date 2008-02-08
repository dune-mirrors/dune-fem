#ifndef DUNE_ISTLMATRIXWRAPPER_HH
#define DUNE_ISTLMATRIXWRAPPER_HH

#if HAVE_DUNE_ISTL

//- system includes 
#include <vector> 

//- Dune istl includes 
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

//- Dune fem includes 
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/function/common/scalarproducts.hh>

namespace Dune { 

  ///////////////////////////////////////////////////////
  // --BlockMatrixHandle
  //////////////////////////////////////////////////////
  template <class LittleBlockType, 
            class RowDiscreteFunctionImp, 
            class ColDiscreteFunctionImp = RowDiscreteFunctionImp> 
  class ImprovedBCRSMatrix : public BCRSMatrix<LittleBlockType> 
  {
    public:
      typedef RowDiscreteFunctionImp RowDiscreteFunctionType;
      typedef ColDiscreteFunctionImp ColDiscreteFunctionType;
      
      typedef BCRSMatrix<LittleBlockType> BaseType; 
      typedef typename BaseType :: RowIterator RowIteratorType ;
      typedef typename BaseType :: ColIterator ColIteratorType ;

      typedef typename BaseType :: size_type size_type;

      //===== type definitions and constants

      //! export the type representing the field
      typedef typename BaseType::field_type field_type;

      //! export the type representing the components
      typedef typename BaseType::block_type block_type;

      //! export the allocator type
      typedef typename BaseType:: allocator_type allocator_type;

      //! implement row_type with compressed vector
      typedef typename BaseType :: row_type row_type;

      //! increment block level counter
      enum {
        //! The number of blocklevels the matrix contains.
        blocklevel = BaseType :: blocklevel 
      };

      /** \brief Iterator for the entries of each row */
      typedef typename BaseType :: ColIterator ColIterator;

      /** \brief Iterator for the entries of each row */
      typedef typename BaseType :: ConstColIterator ConstColIterator;

      /** \brief Const iterator over the matrix rows */
      typedef typename BaseType :: RowIterator RowIterator;

      /** \brief Const iterator over the matrix rows */
      typedef typename BaseType :: ConstRowIterator ConstRowIterator;

      //! type of discrete function space 
      typedef typename ColDiscreteFunctionType :: DiscreteFunctionSpaceType ColSpaceType;

      //! type of row block vector 
      typedef typename RowDiscreteFunctionType :: DofStorageType  RowBlockVectorType; 

      //! type of column block vector 
      typedef typename ColDiscreteFunctionType :: DofStorageType  ColBlockVectorType; 

    protected:  
      size_type nz_;

      int localRows_; 
      int localCols_;

    public:
      //! constructor used by ISTLMatrixObject
      ImprovedBCRSMatrix(size_type rows, size_type cols)
        : BaseType (rows,cols,BaseType::row_wise)
        , nz_(0)
      {
      }

      //! constuctor used by ILU preconditioner 
      ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz)
        : BaseType (rows,cols, BaseType::row_wise)
        , nz_(nz)
      {
      }
      
      //! copy constructor, needed by ISTL preconditioners 
      ImprovedBCRSMatrix(const ImprovedBCRSMatrix& org) 
        : BaseType(org) 
        , nz_(org.nz_)
        , localRows_(org.localRows_)
        , localCols_(org.localCols_)
      {}

      //! setup matrix entires 
      template <class RowMapperType, class ColMapperType,
                class StencilCreatorImp> 
      void setup(const ColSpaceType& colSpace, 
                 const RowMapperType & rowMapper, 
                 const ColMapperType & colMapper,
                 const StencilCreatorImp& stencil, 
                 bool verbose = false) 
      { 
        // if empty grid, do nothing
        if( colSpace.begin() == colSpace.end() ) return ;
        
        {
          // initialize some values 
          localRows_ = rowMapper.maxNumDofs(); 
          localCols_ = colMapper.maxNumDofs(); 

          // map of indices 
          // necessary because element traversal not necessaryly is in
          // ascending order 
          std::map< int , std::set<int> > indices;

          // build matrix entries
          stencil.setup(colSpace, rowMapper, colMapper, indices , (ColDiscreteFunctionType*) 0);

          // type of create interator 
          typedef typename BaseType :: CreateIterator CreateIteratorType; 
          // not insert map of indices into matrix 
          CreateIteratorType endcreate = this->createend();
          for(CreateIteratorType create = this->createbegin();
              create != endcreate; ++create) 
          {
            // set of column indices 
            std::set<int>& localIndices = indices[create.index()];
            typedef typename std::set<int>::iterator iterator;
            iterator end = localIndices.end();
            // insert all indices for this row 
            for (iterator it = localIndices.begin(); it != end; ++it)
            {
              create.insert( *it );
            }
          }
        }

        // in verbose mode some output 
        if(verbose)  
        {
          std::cout << "ISTLMatrix::setup: finished assembly of matrix structure! \n";
        }
      }

      //! clear Matrix, i.e. set all entires to 0
      void clear() 
      {
        RowIteratorType endi=this->end();
        for (RowIteratorType i=this->begin(); i!=endi; ++i)
        {
          ColIteratorType endj = (*i).end();
          for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
          {
            (*j) = 0;
          }
        }
      }

      //! print matrix 
      void print(std::ostream & s) const 
      {
        std::cout << "Print ISTLMatrix \n";
        ConstRowIterator endi=this->end();
        for (ConstRowIterator i=this->begin(); i!=endi; ++i)
        {
          ConstColIterator endj = (*i).end();
          for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          {
            s << (*j) << std::endl;
          }
        }
      }
  };

  template <class RowSpaceImp, class ColSpaceImp, class TraitsImp> 
  class ISTLMatrixObject;
  
  template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
  struct ISTLMatrixTraits
  {
    typedef RowSpaceImp RowSpaceType; 
    typedef ColSpaceImp ColumnSpaceType; 
    typedef ISTLMatrixTraits<RowSpaceType,ColumnSpaceType> ThisType;  
    
    template <class OperatorTraits>
    struct MatrixObject 
    {
      typedef ISTLMatrixObject<RowSpaceType,ColumnSpaceType,OperatorTraits> MatrixObjectType; 
    };
  };

  //! MatrixObject handling an istl matrix 
  template <class RowSpaceImp, class ColSpaceImp, class TraitsImp> 
  class ISTLMatrixObject
  {
  public:  
    //! type of traits 
    typedef TraitsImp Traits;

    //! type of stencil defined by operator
    typedef typename Traits :: StencilType StencilType;
    
    //! type of space defining row structure 
    typedef RowSpaceImp RowSpaceType; 
    //typedef typename Traits :: RowSpaceType RowSpaceType;
    //! type of space defining column structure 
    typedef ColSpaceImp ColumnSpaceType; 
    //typedef typename Traits :: ColumnSpaceType ColumnSpaceType;

    //! type of this pointer 
    typedef ISTLMatrixObject<RowSpaceType,ColumnSpaceType,Traits> ThisType;


  private:  
    typedef typename RowSpaceType::GridType GridType; 
    typedef typename GridType::template Codim<0>::Entity EntityType;

    enum { littleRows = RowSpaceType :: localBlockSize };
    enum { littleCols = ColumnSpaceType :: localBlockSize };
    
    typedef typename RowSpaceType :: RangeFieldType RangeFieldType;
    
    typedef FieldMatrix<RangeFieldType, littleRows, littleCols> LittleBlockType; 

    typedef BlockVectorDiscreteFunction< RowSpaceType >     RowDiscreteFunctionType; 
    typedef typename RowDiscreteFunctionType :: LeakPointerType  RowLeakPointerType;
    typedef BlockVectorDiscreteFunction< ColumnSpaceType >  ColumnDiscreteFunctionType; 
    typedef typename ColumnDiscreteFunctionType :: LeakPointerType  ColumnLeakPointerType;
    
    typedef typename RowDiscreteFunctionType :: DofStorageType    RowBlockVectorType; 
    typedef typename ColumnDiscreteFunctionType :: DofStorageType ColumnBlockVectorType; 

    typedef typename RowSpaceType :: BlockMapperType RowMapperType; 
    typedef typename ColumnSpaceType :: BlockMapperType ColMapperType; 

  public:
    //! type of used matrix 
    typedef ImprovedBCRSMatrix< LittleBlockType , 
                                RowDiscreteFunctionType , 
                                ColumnDiscreteFunctionType > MatrixType;
    typedef typename Traits :: template Adapter < MatrixType > ::  MatrixAdapterType MatrixAdapterType;
    // get preconditioner type from MatrixAdapterType
    typedef ThisType PreconditionMatrixType;
    typedef typename MatrixAdapterType :: ParallelScalarProductType ParallelScalarProductType;
    typedef typename MatrixAdapterType :: CommunicationManagerType CommunicationManagerType;
   
    struct LocalMatrixTraits
    {
      typedef RowSpaceType DomainSpaceType ;
      typedef ColumnSpaceType RangeSpaceType;
      typedef typename RowSpaceType :: RangeFieldType RangeFieldType;
      typedef MatrixType LocalMatrixType;
      typedef typename MatrixType:: block_type LittleBlockType;
    };

    //! LocalMatrix 
    template <class MatrixObjectImp> 
    class LocalMatrix : public LocalMatrixDefault<LocalMatrixTraits>
    {
    public:  
      //! type of base class 
      typedef LocalMatrixDefault<LocalMatrixTraits> BaseType;

      //! type of matrix object 
      typedef MatrixObjectImp MatrixObjectType;
      //! type of matrix 
      typedef typename MatrixObjectImp :: MatrixType MatrixType;
      //! type of little blocks 
      typedef typename MatrixType:: block_type LittleBlockType;
      //! type of entries of little blocks 
      typedef typename RowSpaceType :: RangeFieldType DofType;

      //! type of row mapper 
      typedef typename MatrixObjectType :: RowMapperType RowMapperType;
      //! type of col mapper 
      typedef typename MatrixObjectType :: ColMapperType ColMapperType;
        
    private:
      // special mapper omiting block size 
      const RowMapperType& rowMapper_;
      const ColMapperType& colMapper_;
      
      // number of local matrices 
      int numRows_;
      int numCols_;

      // vector with pointers to local matrices 
      std::vector< std::vector<LittleBlockType *> > matrices_;

      // matrix to build 
      const MatrixObjectType& matrixObj_;

      // type of actual geometry 
      GeometryType geomType_;
      
    public:  
      LocalMatrix(const MatrixObjectType & mObj,
                  const RowSpaceType & rowSpace,
                  const ColumnSpaceType & colSpace)
        : BaseType( rowSpace, colSpace )
        , rowMapper_(mObj.rowMapper())
        , colMapper_(mObj.colMapper())
        , numRows_( rowMapper_.maxNumDofs() )
        , numCols_( colMapper_.maxNumDofs() )
        , matrixObj_(mObj)
        , geomType_(GeometryType::simplex,0)
      {
      }

      void init(const EntityType & rowEntity,
                const EntityType & colEntity)
      {
        if( geomType_ != rowEntity.geometry().type() ) 
        {
          // initialize base functions sets 
          BaseType :: init ( rowEntity , colEntity );

          geomType_ = rowEntity.geometry().type();
          numRows_  = rowMapper_.numDofs(rowEntity);
          numCols_  = colMapper_.numDofs(colEntity);
          matrices_.resize( numRows_ );

          MatrixType& matrix = matrixObj_.matrix();
          for(int i=0; i<numRows_; ++i)
          {
            matrices_[i].resize( numCols_ );
            const int rowIdx = rowMapper_.mapToGlobal(rowEntity,i);
            for(int j=0; j<numCols_; ++j) 
            {
              matrices_[i][j] =
                &matrix[rowIdx][colMapper_.mapToGlobal(colEntity,j)];
            }
          }
        }
        else 
        {
          MatrixType& matrix = matrixObj_.matrix();
          for(int i=0; i<numRows_; ++i)
          {
            const int rowIdx = rowMapper_.mapToGlobal(rowEntity,i);
            for(int j=0; j<numCols_; ++j) 
            {
              matrices_[i][j] =
                &matrix[rowIdx][colMapper_.mapToGlobal(colEntity,j)];
            }
          }
        }
      }

      LocalMatrix(const LocalMatrix& org) 
        : BaseType( org )
        , rowMapper_(org.rowMapper_)
        , colMapper_(org.colMapper_)
        , numRows_( org.numRows_ )
        , numCols_( org.numCols_ ) 
        , matrices_(org.matrices_)
        , matrixObj_(org.matrixObj_)
        , geomType_(org.geomType_)
      {
      }

    private: 
      // check whether given (row,col) pair is valid
      void check(int localRow, int localCol) const 
      {
        const size_t row = (int) localRow / littleRows;
        const size_t col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        assert( row < matrices_.size() ) ;
        assert( col < matrices_[row].size() ); 
        assert( lRow < littleRows );
        assert( lCol < littleCols ); 
      }

      DofType& getValue(const int localRow, const int localCol) 
      {
        const int row = (int) localRow / littleRows;
        const int col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        return (*matrices_[row][col])[lRow][lCol];
      }
    public:
      int rows () const { return matrices_.size()*littleRows; }
      int cols () const { return matrices_[0].size()*littleCols; }

      void add(const int localRow, const int localCol , const DofType value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        getValue(localRow,localCol) += value;
      }

      void set(const int localRow, const int localCol , const DofType value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        getValue(localRow,localCol) = value;
      }

      //! make unit row (all zero, diagonal entry 1.0 )
      void unitRow(const int localRow) 
      {
        const int row = (int) localRow / littleRows;
        const int lRow = localRow%littleRows;

        // get number of columns  
        const int col = cols();
        for(int localCol=0; localCol<col; ++localCol) 
        {
          const int col = (int) localCol / littleCols;
          const int lCol = localCol%littleCols;
          (*matrices_[row][col])[lRow][lCol] = 0;
        }
        // set diagonal entry to 1 
        (*matrices_[row][row])[lRow][lRow] = 1;
      }

      // get entry of matrix 
      DofType get(int localRow, int localCol ) const
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        return getValue(localRow,localCol); 
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        for(int i=0; i<matrices_.size(); ++i)
          for(int j=0; j<matrices_[i].size(); ++j)
            (*matrices_[i][j]) = (DofType) 0;
      }

      //! empty as the little matrices are already sorted
      void resort() 
      {
      }
    };

  public:
    //! type of local matrix 
    typedef LocalMatrix<ThisType> ObjectType;
    typedef ThisType LocalMatrixFactoryType;
    typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
    //! type of local matrix 
    typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

  private:  
    const RowSpaceType & rowSpace_;
    const ColumnSpaceType & colSpace_;

    // sepcial row mapper 
    RowMapperType& rowMapper_;
    // special col mapper 
    ColMapperType& colMapper_;

    int size_;

    int sequence_;

    mutable MatrixType* matrix_;

    CommunicationManagerType comm_;
    ParallelScalarProductType scp_;

    int numIterations_; 
    double relaxFactor_; 
      
    enum PreConder_Id { none  = 0 , // no preconditioner 
                        ssor  = 1 , // SSOR preconditioner 
                        sor   = 2 , // SOR preconditioner 
                        ilu_0 = 3 , // ILU-0 preconditioner 
                        ilu_n = 4 , // ILU-n preconditioner 
                        gauss_seidel= 5 , // Gauss-Seidel preconditioner 
                        jacobi = 6  // Jacobi preconditioner 
    };
    
    PreConder_Id preconditioning_;

    mutable LocalMatrixStackType localMatrixStack_;

    mutable MatrixAdapterType* matrixAdap_;
    mutable RowBlockVectorType* Arg_;
    mutable ColumnBlockVectorType* Dest_;

    // prohibit copy constructor 
    ISTLMatrixObject(const ISTLMatrixObject&); 
  public:  
    //! constructor 
    //! \param rowSpace space defining row structure 
    //! \param colSpace space defining column structure 
    //! \param paramfile parameter file to read variables 
    //!         - Preconditioning: {0,1,2,3,4,5,6} put -1 to get info
    //!         - Pre-iteration: number of iteration of preconditioner
    //!         - Pre-relaxation: relaxation factor   
    ISTLMatrixObject ( const RowSpaceType &rowSpace,
                       const ColumnSpaceType &colSpace,
                       const std :: string &paramfile = "" )
      : rowSpace_(rowSpace)
      , colSpace_(colSpace)
      // create scp to have at least one instance 
      // otherwise instance will be deleted during setup
      // get new mappers with number of dofs without considerung block size 
      , rowMapper_( rowSpace.blockMapper() )
      , colMapper_( colSpace.blockMapper() )
      , size_(-1)
      , sequence_(-1)
      , matrix_(0)
      , comm_(rowSpace_)
      , scp_(colSpace_)
      , numIterations_(5)
      , relaxFactor_(1.1)
      , preconditioning_(none)
      , localMatrixStack_( *this )
      , matrixAdap_(0)
      , Arg_(0)
      , Dest_(0)
    {
      if(paramfile != "")
      {
        const bool output = (rowSpace_.grid().comm().rank() == 0);
        int preCon = 0;
        readParameter(paramfile,"Preconditioning",preCon, output);
        if( preCon >= 0 && preCon <= 6) 
          preconditioning_ = (PreConder_Id) preCon;
        else 
          preConErrorMsg(preCon);
        
        readParameter(paramfile,"Pre-iteration",numIterations_, output);
        readParameter(paramfile,"Pre-relaxation",relaxFactor_, output);
      }

      assert( rowMapper_.size() == colMapper_.size() );
    }

  public:  
    //! destructor 
    ~ISTLMatrixObject() 
    {
      removeObj();
    }

    //! return reference to system matrix 
    MatrixType & matrix() const 
    { 
      assert( matrix_ );
      return *matrix_; 
    }
    
    //! return matrix adapter object  
    MatrixAdapterType matrixAdapter() const 
    { 
      // no preconditioner 
      if( preconditioning_ == none )
      {
        return MatrixAdapterType(matrix(),rowSpace_,colSpace_);
      }
      // SSOR 
      else if( preconditioning_ == ssor )
      {
        typedef SeqSSOR<MatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return MatrixAdapterType(matrix(),rowSpace_,colSpace_,
                                 numIterations_ , relaxFactor_,(PreconditionerType*)0);
      }
      // SOR 
      else if(preconditioning_ == sor )
      {
        typedef SeqSOR<MatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return MatrixAdapterType(matrix(),rowSpace_,colSpace_,
                                 numIterations_ , relaxFactor_,(PreconditionerType*)0);
      }
      // ILU-0 
      else if(preconditioning_ == ilu_0)
      {
        typedef SeqILU0<MatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return MatrixAdapterType(matrix(),rowSpace_,colSpace_,
                                 relaxFactor_,(PreconditionerType*)0);
      }
      // ILU-n
      else if(preconditioning_ == ilu_n)
      {
        typedef SeqILUn<MatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return MatrixAdapterType(matrix(),rowSpace_,colSpace_,
                                 numIterations_ , relaxFactor_,(PreconditionerType*)0);
      }
      // Gauss-Seidel
      else if(preconditioning_ == gauss_seidel)
      {
        typedef SeqGS<MatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return MatrixAdapterType(matrix(),rowSpace_,colSpace_,
                                 numIterations_ , relaxFactor_,(PreconditionerType*)0);
      }
      // Jacobi 
      else if(preconditioning_ == jacobi)
      {
        typedef SeqJac<MatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return MatrixAdapterType(matrix(),rowSpace_,colSpace_,
                                 numIterations_ , relaxFactor_,(PreconditionerType*)0);
      }
      else 
      {
        preConErrorMsg(preconditioning_);
      }
      return MatrixAdapterType(matrix(),rowSpace_,colSpace_ );
    }
    
    //! return true, because in case of no preconditioning we have empty
    //! preconditioner (used by OEM methods)
    bool hasPreconditionMatrix() const { return (preconditioning_ != none); }

    //! return reference to preconditioner object (used by OEM methods)
    const PreconditionMatrixType& preconditionMatrix() const { return *this; }

    //! set all matrix entries to zero 
    void clear()
    {
      matrix().clear();
    }

    //! reserve matrix with right size 
    void reserve(bool verbose = false) 
    {
      // if grid sequence number changed, rebuild matrix 
      if(sequence_ != rowSpace_.sequence())
      {
        removeObj();

        matrix_ = new MatrixType(rowMapper_.size(), colMapper_.size());
        StencilType stencil; 
        matrix().setup(colSpace_,rowMapper(),colMapper(),stencil,verbose);

        sequence_ = rowSpace_.sequence();
      }
    }

    //! we only have right precondition
    bool rightPrecondition() const { return true; }

    //! precondition method for OEM Solvers 
    //! not fast but works, double is copied to block vector 
    //! and after application copied back
    void precondition(const double* arg, double* dest) const
    {
      createBlockVectors();
      
      assert( Arg_ );
      assert( Dest_ );

      RowBlockVectorType& Arg = *Arg_;
      ColumnBlockVectorType & Dest = *Dest_;
      
      // copy from double 
      double2Block(arg, Arg);
      
      // set Dest to zero 
      Dest = 0;
      
      assert( matrixAdap_ );
      // not parameter swaped for preconditioner 
      matrixAdap_->preconditionAdapter().apply(Dest , Arg);

      // copy back 
      block2Double( Dest , dest);
    }
      
    //! mult method for OEM Solver 
    void multOEM(const double* arg, double* dest) const
    {
      createBlockVectors();
      
      assert( Arg_ );
      assert( Dest_ );

      RowBlockVectorType& Arg = *Arg_;
      ColumnBlockVectorType & Dest = *Dest_;
      
      // copy from double 
      double2Block(arg, Arg);
      
      // call mult of matrix adapter  
      assert( matrixAdap_ );
      matrixAdap_->apply( Arg, Dest );

      //  copy back 
      block2Double( Dest , dest);
    }
    
    //! apply with discrete functions 
    void apply(const RowDiscreteFunctionType& arg,
               RowDiscreteFunctionType& dest) const 
    {
      createMatrixAdapter();
      assert( matrixAdap_ );
      matrixAdap_->apply( arg.blockVector(), dest.blockVector() );
    }

    //! apply with arbitrary discrete functions calls multOEM 
    template <class RowDFType, class ColDFType>
    void apply(const RowDFType& arg, ColDFType& dest) const 
    {
      multOEM( arg.leakPointer(), dest.leakPointer ());
    }

    //! mult method of matrix object used by oem solver
    void multOEM(const RowLeakPointerType& arg, ColumnLeakPointerType& dest) const
    {
      createMatrixAdapter();
      assert( matrixAdap_ );
      matrixAdap_->apply( arg.blockVector(), dest.blockVector() );
    }

    //! resort row numbering in matrix to have ascending numbering 
    void resort()
    {
    }

    //! create precondition matrix does nothing because preconditioner is
    //! created only when requested 
    void createPreconditionMatrix() 
    {
    }

    //! print matrix 
    void print(std::ostream & s) const 
    { 
      matrix().print(std::cout);
    }

    const RowMapperType& rowMapper() const { return rowMapper_; }
    const ColMapperType& colMapper() const { return colMapper_; }

    //! interface method from LocalMatrixFactory 
    ObjectType* newObject() const 
    {
      return new ObjectType(*this,
                            rowSpace_,
                            colSpace_);
    }

    //! return local matrix object 
    LocalMatrixType localMatrix(const EntityType& rowEntity, 
                                const EntityType& colEntity) const 
    {
      return LocalMatrixType(localMatrixStack_,rowEntity,colEntity);
    }

  protected:  
    void preConErrorMsg(int preCon) const 
    {
      std::cerr << "ERROR: Wrong precoditioning number (p = " << preCon;
      std::cerr <<") in ISTLMatrixObject! \n";
      std::cerr <<"Valid values are: \n";
      std::cerr <<"0 == no \n";
      std::cerr <<"1 == SSOR \n";
      std::cerr <<"2 == SOR \n";
      std::cerr <<"3 == ILU-0 \n";
      std::cerr <<"4 == ILU-n \n";
      std::cerr <<"5 == Gauss-Seidel \n";
      std::cerr <<"6 == Jacobi \n";
      assert(false);
      exit(1);
    }

    void removeObj() 
    {
      delete Dest_; Dest_ = 0;
      delete Arg_;  Arg_ = 0;
      delete matrixAdap_; matrixAdap_ = 0;
      delete matrix_; matrix_ = 0;
    }

    // copy double to block vector 
    void double2Block(const double* arg, RowBlockVectorType& dest) const 
    {
      typedef typename RowBlockVectorType :: block_type BlockType;
      const size_t blocks = dest.size();
      int idx = 0;
      for(size_t i=0; i<blocks; ++i) 
      {
        BlockType& block = dest[i];
        enum { blockSize = BlockType :: dimension };
        for(int j=0; j<blockSize; ++j, ++idx) 
        {
          block[j] = arg[idx];
        }
      }
    }
      
    // copy block vector to double 
    void block2Double(const ColumnBlockVectorType& arg, double* dest) const 
    {
      typedef typename ColumnBlockVectorType :: block_type BlockType;
      const size_t blocks = arg.size();
      int idx = 0;
      for(size_t i=0; i<blocks; ++i) 
      {
        const BlockType& block = arg[i];
        enum { blockSize = BlockType :: dimension };
        for(int j=0; j<blockSize; ++j, ++idx) 
        {
          dest[idx] = block[j];
        }
      }
    }

    void createBlockVectors() const
    {
      if( ! Arg_ || ! Dest_ ) 
      { 
        delete Arg_; delete Dest_;
        Arg_  = new RowBlockVectorType( rowMapper_.size() );
        Dest_ = new ColumnBlockVectorType( colMapper_.size() );
      }

      createMatrixAdapter ();
    }

    void createMatrixAdapter () const 
    {
      if( ! matrixAdap_ ) 
      { 
        matrixAdap_ = new MatrixAdapterType(matrixAdapter());
      }
    }
    
  };

} // end namespace Dune 
#endif

#endif
