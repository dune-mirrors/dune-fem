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

namespace Dune { 

  ///////////////////////////////////////////////////////
  // --BlockMatrixHandle
  //////////////////////////////////////////////////////
  template <class LittleBlockType, class DiscreteFunctionType> 
  class ImprovedBCRSMatrix : public BCRSMatrix<LittleBlockType> 
  {
    public:
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
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType SpaceType;

      //! type of block vector 
      typedef typename DiscreteFunctionType :: DofStorageType  BlockVectorType; 

      //! type of used communication manager  
      typedef CommunicationManager<SpaceType> CommunicationManagerType; 

    private:  
      size_type nz_;

      int localRows_; 
      int localCols_;

      //! our function space, needed for communication  
      const SpaceType* space_; 

      //! communication manager 
      mutable CommunicationManagerType* comm_;

      std::vector<int> overlapRows_;

    public:
      //! constructor used by ISTLMatrixObject
      ImprovedBCRSMatrix(const SpaceType & space, 
                         CommunicationManagerType& comm,
                         size_type rows, size_type cols)
        : BaseType (rows,cols,BaseType::row_wise)
        , nz_(0)
        , space_(&space)
        , comm_(&comm)         
      {
        // only works for non-hybrid grids so far 
        if( space_->multipleGeometryTypes() )
        {
          DUNE_THROW(NotImplemented,"ISTLMatrix::setup: use of matrix for hybrid grids not implemented yet!");
        }

        //std::cout << "Create Matrix with " << rows << " x " << cols << "\n";
      }

      //! constuctor used by ILU preconditioner 
      ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz)
        : BaseType (rows,cols, BaseType::row_wise)
        , nz_(nz)
        , space_(0)
        , comm_(0)         
      {
        //std::cout << "Create Matrix with " << rows << " x " << cols << "\n";
      }
      
      //! copy constructor, needed by ISTL preconditioners 
      ImprovedBCRSMatrix(const ImprovedBCRSMatrix& org) 
        : BaseType(org) 
        , nz_(org.nz_)
        , localRows_(org.localRows_)
        , localCols_(org.localCols_)
        , space_(org.space_)
        , comm_(org.comm_)
        , overlapRows_(org.overlapRows_) 
      {}

      //! matrix multiplication for OEM solvers 
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

      //! setup matrix entires 
      template <class GridPartType, 
                class RowMapperType, class ColMapperType,
                class StencilCreatorImp> 
      void setup(const GridPartType& gridPart, 
                 const RowMapperType & rowMapper, 
                 const ColMapperType & colMapper,
                 const StencilCreatorImp& stencil, 
                 bool verbose = false) 
      { 
        // get size estimate   
        int size = rowMapper.size();
        size = (int) size / 10;
        overlapRows_.reserve( size );
        overlapRows_.resize(0);

        // if empty grid, do nothing
        if( gridPart.template begin<0> () == gridPart.template end<0> () ) return ;
        {
          // initialize some values 
          localRows_ = rowMapper.numDofs(); 
          localCols_ = colMapper.numDofs(); 

          // map of indices 
          // necessary because element traversal not necessaryly is in
          // ascending order 
          std::map< int , std::set<int> > indices;

          // call stencil creator 
          stencil.setup(gridPart,rowMapper,colMapper,indices,overlapRows_);

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

        // sort overlap rows because of cache efficiency 
        std::sort(overlapRows_.begin(), overlapRows_.end());

        // in verbose mode some output 
        if(verbose)  
        {
          std::cout << "ISTLMatrix::setup: finished assembly of matrix structure! \n";
        }
      }

      //! clear Matrix, i.e. set all entires to 0
      void clear() 
      {
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

        if( LittleBlockType :: rows == LittleBlockType :: cols )
        {
          // for non-interior entities set diag to 1 for ILU Preconditioner 
          const int overlap = overlapRows_.size();
          for(int i=0; i<overlap; ++i)
          {
            const int idx = overlapRows_[i]; 
            LittleBlockType& diag = this->operator[](idx)[idx];
            for(int k=0; k<LittleBlockType :: rows; ++k)  
            {
              diag[k][k] = 1;
            }
          }
        }
      }

      //! print matrix 
      void print(std::ostream & s) const 
      {
        std::cout << "Print ISTLMatrix \n";
        RowIteratorType endi=this->end();
        for (RowIteratorType i=this->begin(); i!=endi; ++i)
        {
          ColIteratorType endj = (*i).end();
          for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
          {
            s << (*j) << std::endl;
          }
        }
      }

      //! communicate block vector 
      void communicate(const BlockVectorType& arg) const 
      {
        if(comm_)
        {
          assert( space_ );
          // if serial run, just return 
          if(space_->grid().comm().size() <= 1) 
          {
            return;
          }

          // exchange data 
          DiscreteFunctionType tmp("ImprovedBCRSMatrix::communicate_tmp",*space_,arg);
          comm_->exchange( tmp );
        }
      }

      //! apply matrix: \f$ y = A(x) \f$
      void mult(const BlockVectorType& x, BlockVectorType& y) const 
      {
        // exchange data 
        communicate( x );

        // clear vector  
        y = 0;
        // multiply 
        this->umv(x,y);

        // delete non interior entries 
        deleteNonInterior(y);
      }

      //! apply scaled: \f$ y = y + \alpha A(x) \f$
      void multAdd(field_type alpha, const BlockVectorType& x, BlockVectorType& y) const 
      {
        // exchange data 
        communicate( x );

        this->usmv(alpha,x,y);

        // delete non interior entries 
        deleteNonInterior(y);
      }

    private:  
      // delete all vector entries that belong not to interior entities 
      void deleteNonInterior(BlockVectorType& y) const
      {
        // set all entries belonging to non-interior elements to zero 
        const int overlap = overlapRows_.size();
        for(int i=0; i<overlap; ++i)
        {
          y[overlapRows_[i]] = 0;
        }
      }
  };

  //! wrapper class to store perconditioner 
  //! as the interface class does not have to category 
  //! enum 
  template<class X, class Y>
  class PreconditionerWrapper : public Preconditioner<X,Y>
  {
    typedef Preconditioner<X,Y> PreconditionerInterfaceType;
    PreconditionerInterfaceType* preconder_; 
    
    //! set preconder to zero 
    PreconditionerWrapper (const PreconditionerWrapper&);
  public:
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    enum {
      //! \brief The category the precondtioner is part of.
      category=SolverCategory::sequential};

    //! set preconder to zero 
    PreconditionerWrapper () : preconder_(0) {}
    
    //! create preconditioner of given type 
    template <class MatrixType, class PreconditionerType>
    PreconditionerWrapper(MatrixType & m, int iter, field_type relax, const PreconditionerType*) 
      : preconder_(new PreconditionerType(m,iter,relax)) {}
    
    //! create preconditioner of given type 
    template <class MatrixType, class PreconditionerType>
    PreconditionerWrapper(MatrixType & m, field_type relax, const PreconditionerType*) 
      : preconder_(new PreconditionerType(m,relax)) {}
    
    //! \copydoc Preconditioner 
    virtual void pre (X& x, Y& b) 
    {
      // apply preconditioner
      if( preconder_ ) preconder_->pre(x,b);
    }

    //! \copydoc Preconditioner 
    virtual void apply (X& v, const Y& d)
    {
      if( preconder_ ) 
      {
        // apply preconditioner
        preconder_->apply(v,d);
      }
      else 
      {
        // just copy values 
        v = d;
      }
    }

    //! \copydoc Preconditioner 
    virtual void post (X& x) {
      // apply preconditioner
      if( preconder_ ) preconder_->post(x);
    }

    // every abstract base class has a virtual destructor
    virtual ~PreconditionerWrapper () 
    {
      delete preconder_;
      preconder_ = 0;
    }
  };

  //! MatrixObject handling an istl matrix 
  template <class RowSpaceImp, class ColumnSpaceImp> 
  class ISTLMatrixObject  
  {
  public:  
    //! type of space defining row structure 
    typedef RowSpaceImp RowSpaceType;
    //! type of space defining column structure 
    typedef ColumnSpaceImp ColumnSpaceType;

    //! type of this pointer 
    typedef ISTLMatrixObject<RowSpaceImp,ColumnSpaceImp> ThisType;

  private:  
    typedef typename RowSpaceType::GridType GridType; 
    typedef typename GridType::template Codim<0>::Entity EntityType;

    enum { littleRows = RowSpaceType :: localBlockSize };
    enum { littleCols = ColumnSpaceType :: localBlockSize };
    
    typedef typename RowSpaceType :: RangeFieldType RangeFieldType;
    
    typedef FieldMatrix<RangeFieldType, littleRows, littleCols> LittleBlockType; 

    typedef BlockVectorDiscreteFunction< RowSpaceType >  DiscreteFunctionType; 
    typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType; 

    typedef typename RowSpaceType :: BlockMapperType RowMapperType; 
    typedef typename ColumnSpaceType :: BlockMapperType ColMapperType; 

  public:
    //! type of used matrix 
    typedef ImprovedBCRSMatrix< LittleBlockType , DiscreteFunctionType > MatrixType;
   
    //! type of preconditioner 
    typedef PreconditionerWrapper<BlockVectorType,BlockVectorType> PreconditionMatrixType;

    struct LocalMatrixTraits
    {
      typedef RowSpaceImp DomainSpaceType ;
      typedef ColumnSpaceImp RangeSpaceType;
      typedef typename RowSpaceImp :: RangeFieldType RangeFieldType;
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
        , numRows_( rowMapper_.numDofs() )
        , numCols_( colMapper_.numDofs() )
        , matrixObj_(mObj)
        , geomType_(GeometryType::simplex,0)
      {
      }

      void init(const EntityType & rowEntity,
                const EntityType & colEntity)
      {
        if( geomType_ != rowEntity.geometry().type() ) 
        {
          geomType_ = rowEntity.geometry().type();
          numRows_  = rowMapper_.numDofs();
          numCols_  = colMapper_.numDofs();
          matrices_.resize( numRows_ );

          MatrixType& matrix = matrixObj_.matrix();
          //std::cout << matrices_.size() << " size rows \n";
          for(int i=0; i<numRows_; ++i)
          {
            matrices_[i].resize( numCols_ );
            //std::cout << matrices_[i].size() << " size col \n";
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
        return 0.0;
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        return getValue(localRow,localCol); 
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        for(int i=0; i<matrices_.size(); ++i)
        {
          (*matrices_[i]) = (DofType) 0;
        }
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
    typedef CommunicationManager<RowSpaceType> CommunicationManagerType;

    const RowSpaceType & rowSpace_;
    const ColumnSpaceType & colSpace_;

    // sepcial row mapper 
    RowMapperType& rowMapper_;
    // special col mapper 
    ColMapperType& colMapper_;

    int size_;

    int sequence_;

    mutable MatrixType* matrix_;
    mutable PreconditionMatrixType* preconder_;

    CommunicationManagerType comm_;

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
    ISTLMatrixObject(const RowSpaceType & rowSpace,
                     const ColumnSpaceType & colSpace,
                     const std::string& paramfile)
      : rowSpace_(rowSpace)
      , colSpace_(colSpace)
      // get new mappers with number of dofs without considerung block size 
      , rowMapper_( rowSpace.blockMapper() )
      , colMapper_( colSpace.blockMapper())
      , size_(-1)
      , sequence_(-1)
      , matrix_(0)
      , preconder_(0)
      , comm_(rowSpace_)
      , numIterations_(5)
      , relaxFactor_(1.1)
      , preconditioning_(none)
      , localMatrixStack_( *this )
    {
      if(paramfile != "")
      {
        int preCon = 0;
        readParameter(paramfile,"Preconditioning",preCon);
        if( preCon >= 0 && preCon <= 6) 
          preconditioning_ = (PreConder_Id) preCon;
        else 
          preConErrorMsg(preCon);
        
        readParameter(paramfile,"Pre-iteration",numIterations_);
        readParameter(paramfile,"Pre-relaxation",relaxFactor_);
      }

      // only ILU-0 works savely in parallel
      if(rowSpace_.grid().comm().size() > 1)
      {
        if( (preconditioning_ != none) && (preconditioning_ != ilu_0) )
        {
          std::cerr << "ERROR: Only Preconditioner ILU-0 works in parallel! " << std::endl;
          abort();
        }
      }
      
      assert( rowMapper_.size() == colMapper_.size() );
    }

    //! destructor 
    ~ISTLMatrixObject() 
    {
      delete preconder_;
      delete matrix_;
    }

    //! return reference to system matrix 
    MatrixType & matrix() const 
    { 
      assert( matrix_ );
      return *matrix_; 
    }
    
    //! return true, because in case of no preconditioning we have empty
    //! preconditioner 
    bool hasPcMatrix () const { return true; }

    //! return reference to preconditioner
    const PreconditionMatrixType& pcMatrix () const  
    { 
      if( !preconder_ )
      {
        preconder_ = createPreconditioner();
      }
      return *preconder_; 
    }

    //! set all matrix entries to zero 
    void clear()
    {
      matrix().clear();
    }

    //! reserve matrix with right size 
    template <class StencilCreatorImp>
    void reserve(const StencilCreatorImp& stencil, 
                 bool verbose = false) 
    {
      // if grid sequence number changed, rebuild matrix 
      if(sequence_ != rowSpace_.sequence())
      {
        delete matrix_; matrix_ = 0;
        delete preconder_; preconder_ = 0;

        matrix_ = new MatrixType(rowSpace_, comm_, rowMapper_.size(), colMapper_.size());
        matrix().setup(rowSpace_.gridPart(),rowMapper(),colMapper(),stencil,verbose);

        sequence_ = rowSpace_.sequence();
      }
    }

    //! mult method of matrix object used by oem solver
    void multOEM(const double * arg, double * dest) const
    {
      matrix().multOEM(arg,dest);
    }

    //! mult method of matrix object used by oem solver
    template <class LeakPtrImp> 
    void multOEM(const LeakPtrImp& arg, LeakPtrImp& dest) const
    {
      DUNE_THROW(NotImplemented,"ISTLMatrixObject::multOEM for arbitrary type not implemented!");
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

    LocalMatrixType localMatrix(const EntityType& rowEntity, 
                                const EntityType& colEntity) const 
    {
      return LocalMatrixType(localMatrixStack_,rowEntity,colEntity);
    }

  private:  
    PreconditionMatrixType* createPreconditioner() const
    {
      // no preconditioner 
      if( preconditioning_ == none )
      {
        return new PreconditionMatrixType();
      }
      // SSOR 
      else if( preconditioning_ == ssor )
      {
        typedef SeqSSOR<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // SOR 
      else if(preconditioning_ == sor )
      {
        typedef SeqSOR<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // ILU-0 
      else if(preconditioning_ == ilu_0)
      {
        typedef SeqILU0<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), relaxFactor_, (PreconditionerType*)0);
      }
      // ILU-n
      else if(preconditioning_ == ilu_n)
      {
        typedef SeqILUn<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // Gauss-Seidel
      else if(preconditioning_ == gauss_seidel)
      {
        typedef SeqGS<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // Jacobi 
      else if(preconditioning_ == jacobi)
      {
        typedef SeqJac<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      else 
      {
        preConErrorMsg(preconditioning_);
      }
      return 0;
    }

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
  };

} // end namespace Dune 
#endif

#endif
