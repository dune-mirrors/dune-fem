#ifndef DUNE_ISTLMATRIXWRAPPER_HH
#define DUNE_ISTLMATRIXWRAPPER_HH

//- system includes 
#include <vector> 

//- Dune istl includes 
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

//- Dune fem includes 
#include <dune/fem/discretefunction/staticfunction.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune { 

//#if HAVE_DUNE_ISTL 
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

      size_type nz_;

      int localRows_; 
      int localCols_;

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

      //! our function space  
      const SpaceType* space_; 

      //! communication manager 
      mutable CommunicationManagerType* comm_;

      std::vector<int> overlapRows_;
    public:
      //! constructor used by ISTLMatrixObject
      ImprovedBCRSMatrix(const SpaceType & space, 
                         CommunicationManagerType& comm,
                         size_type rows, size_type cols, size_type nz)
        : BaseType (rows,cols, BaseType::row_wise)
        , nz_(nz)
        , space_(&space)
        , comm_(&comm)         
      {
        std::cout << "Create Matrix with " << rows << " x " << cols << "\n";
      }

      //! constuctor used by ILU preconditioner 
      ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz)
        : BaseType (rows,cols, BaseType::row_wise)
        , nz_(nz)
        , space_(0)
        , comm_(0)         
      {
        std::cout << "Create Matrix with " << rows << " x " << cols << "\n";
      }

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
      template <class RowSpaceType, class ColSpaceType,
                class GridType,
                PartitionIteratorType pitype,
                template <class,PartitionIteratorType> class GridPartType> 
      void setup(const RowSpaceType & rowSpace, 
                 const ColSpaceType & colSpace,
                 const GridPartType<GridType,pitype> & gridPart) 
      { 
        int size = rowSpace.indexSet().size(0);
        size = (int) size / 10;
        overlapRows_.reserve( size );
        overlapRows_.resize(0);
        {
          
          typedef typename BaseType :: CreateIterator CreateIteratorType; 

          CreateIteratorType create = this->createbegin();
          CreateIteratorType endcreate = this->createend();

          //! we need all partition iterator here  
          typedef GridPartType<GridType,All_Partition> AllPartType; 
          typedef typename AllPartType :: template Codim<0> :: IteratorType  IteratorType;
          AllPartType allPart(const_cast<GridType&> (rowSpace.grid()));
          
          typedef typename GridType :: template Codim<0> :: Entity EntityType;
          typedef typename AllPartType:: IntersectionIteratorType IntersectionIteratorType; 
          
          IteratorType endit = allPart.template end<0> ();
          for(IteratorType it = allPart.template begin<0> (); it != endit; ++it)
          {
            assert( create != endcreate );

            EntityType & en = *it;

            localRows_ = rowSpace.getBaseFunctionSet(en).numBaseFunctions();
            localCols_ = colSpace.getBaseFunctionSet(en).numBaseFunctions();

            const int elIndex = rowSpace.indexSet().index(en);
            create.insert( elIndex );
            
            if(en.partitionType() != InteriorEntity) 
              overlapRows_.push_back( elIndex );

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
        clear();
        std::sort(overlapRows_.begin(), overlapRows_.end());
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
              (*j) = 0.0;
            }
          }
        }

        // for non-interior entities set diag to 1 for ILU Preconditioner 
        const int overlap = overlapRows_.size();
        for(int i=0; i<overlap; ++i)
        {
          const int idx = overlapRows_[i]; 
          LittleBlockType& diag = this->operator[](idx)[idx];
          CompileTimeChecker<LittleBlockType :: rows == LittleBlockType :: cols > ();
          for(int k=0; k<LittleBlockType :: rows; ++k)  
          {
            diag[k][k] = 1.0;
          }
        }
      }

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

        // multiply 
        y = 0;
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

  template <class RowSpaceType, class ColumnSpaceType> 
  class ISTLMatrixObject  
  {
    typedef typename RowSpaceType::GridType GridType; 
    typedef typename GridType::template Codim<0>::Entity EntityType;
  public:  
    enum { littleRows = RowSpaceType :: numBaseFunctions };
    enum { littleCols = ColumnSpaceType :: numBaseFunctions };
    
    typedef typename RowSpaceType :: RangeFieldType RangeFieldType;
    
    typedef FieldMatrix<RangeFieldType, littleRows, littleCols> LittleBlockType; 
    typedef BlockVector< FieldVector<RangeFieldType, littleRows> > BlockVectorType; 

    typedef StaticDiscreteFunction< RowSpaceType , BlockVectorType>  DiscreteFunctionType; 
    typedef ImprovedBCRSMatrix< LittleBlockType , DiscreteFunctionType > MatrixType;
   
    typedef SeqILUn<MatrixType,BlockVectorType,BlockVectorType> PreconditionMatrixType;

    template <class MatrixImp> 
    class LocalMatrix
    {
      typedef MatrixImp MatrixType;
      typedef typename MatrixType:: block_type LittleBlockType;
      
      const int rowIndex_;
      const int colIndex_;

      LittleBlockType & matrix_;
      
    public:  
      LocalMatrix(MatrixType & m,
                  const EntityType & rowEntity,
                  const RowSpaceType & rowSpace,
                  const EntityType & colEntity,
                  const ColumnSpaceType & colSpace)
        : rowIndex_(rowSpace.indexSet().index(rowEntity))
        , colIndex_(colSpace.indexSet().index(colEntity))
        , matrix_(m[rowIndex_][colIndex_])
      {
      }

    private: 
      LocalMatrix(const LocalMatrix&);

    public:
      void add(int localRow, int localCol , const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );
        matrix_[localRow][localCol] += value;
      }

      void set(int localRow, int localCol , const double value)
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );
        matrix_[localRow][localCol] = value;
      }

      double get(int localRow, int localCol ) const
      {
        assert( localRow >= 0 );
        assert( localCol >= 0 );
        return matrix_[localRow][localCol];
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        matrix_ = 0.0;
      }

      void resort() 
      {
      }
    };

    //! default number of non-zeros per row 
    template <class GridImp>
    struct NonZeros
    {
      enum { dimension = GridImp::dimension };
      static int bound() 
      {
        // upper bound is number of max number of neighbors + myself
        return (2 * dimension) + 1;
      }
    };

    //! default number of non-zeros per row 
    template <int dim, int dimworld>
    struct NonZeros<ALUSimplexGrid<dim,dimworld> > 
    {
      enum { dimension = dim };
      static int bound() 
      {
        // we have dim+1 neighbors when conforming
        // we have 2 * dim-1 (4 or 2) neihbors on each face when
        // non-conforming
        // upper bound is number of max number of neighbors + myself
        return (dimension+1) * (2 * dimension-1) + 1;
      }
    };

    //! default number of non-zeros per row 
    template <int dim, int dimworld>
    struct NonZeros<ALUCubeGrid<dim,dimworld> >
    {
      enum { dimension = dim };
      static int bound() 
      {
        // we have 2 * dim neighbors when conforming
        // we have 2 * dim-1 (4 or 2) neihbors on each face when
        // non-conforming
        // upper bound is number of max number of neighbors + myself
        return (2* dimension) * (2 * dimension-1) + 1;
      }
    };

  public:
    typedef LocalMatrix<MatrixType> LocalMatrixType;
    typedef CommunicationManager<RowSpaceType> CommunicationManagerType;

    const RowSpaceType & rowSpace_;
    const ColumnSpaceType & colSpace_;
    const int factor_; 

    int size_;

    mutable MatrixType* matrix_;
    mutable PreconditionMatrixType* preconder_;

    CommunicationManagerType comm_;

    int numIterations_; 
    double relaxFactor_; 
      
    const bool preconditioning_;

    //! setup matrix handler 
    ISTLMatrixObject(const RowSpaceType & rowSpace,
                     const ColumnSpaceType & colSpace,
                     bool preconditioning,
                     const std::string& paramfile)
      : rowSpace_(rowSpace)
      , colSpace_(colSpace)
      , factor_( NonZeros<GridType>::bound() )
      , size_(-1)
      , matrix_(0)
      , preconder_(0)
      , comm_(rowSpace_)
      , numIterations_(5)
      , relaxFactor_(1.1)
      , preconditioning_(preconditioning)
    {
      if(paramfile != "")
      {
        readParameter(paramfile,"ILU-iteration",numIterations_);
        readParameter(paramfile,"ILU-relaxation",relaxFactor_);
      }
      
      assert( rowSpace_.indexSet().size(0) ==
              colSpace_.indexSet().size(0) );
      reserve(true);
    }

    MatrixType & matrix() const 
    { 
      assert( matrix_ );
      return *matrix_; 
    }
    
    //! return true if precoditioning matrix is provided 
    bool hasPcMatrix () const { return preconditioning_; }
    const PreconditionMatrixType& pcMatrix () const  
    { 
      if( !preconder_ )
      {
        preconder_ = new PreconditionMatrixType( matrix() , 
                             numIterations_ , relaxFactor_ );
      }
      return *preconder_; 
    }

    bool hasBeenSetup () const
    {
      return (size_ > 0);
    }

    void clear()
    {
      matrix().clear();
    }

    void resize(bool verbose = false)
    {
      reserve(verbose);
    }

    void reserve(bool verbose = false) 
    {
      {
        delete matrix_;
        delete preconder_; preconder_ = 0;
        size_ = rowSpace_.indexSet().size(0);

        /*
        { 
          typedef typename RowSpaceType :: IteratorType IteratorType;
          int rows = 0;
          IteratorType endit = rowSpace_.end();
          for(IteratorType it = rowSpace_.begin(); it != endit; ++it)
          {
            ++rows;   
          }
          size_ = rows;
        }
        */
        
        matrix_ = new MatrixType(rowSpace_, comm_, size_, colSpace_.indexSet().size(0),factor_);
        matrix().setup(rowSpace_,colSpace_,rowSpace_.gridPart());
      }
    }

    //! mult method of matrix object used by oem solver
    void multOEM(const double * arg, double * dest) const
    {
    }

    //! resort row numbering in matrix to have ascending numbering 
    void resort()
    {
    }

    void createPreconditionMatrix() 
    {
    }

    void print(std::ostream & s) const 
    { 
      matrix().print(std::cout);
    }
  };
//#endif // HAVE_DUNE_ISTL

} // end namespace Dune 
#endif
