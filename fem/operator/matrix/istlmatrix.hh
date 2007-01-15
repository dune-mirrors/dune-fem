#ifndef DUNE_ISTLMATRIXWRAPPER_HH
#define DUNE_ISTLMATRIXWRAPPER_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/fem/discretefunction/staticfunction.hh>

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
    public:
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

      ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz)
        : BaseType (rows,cols, BaseType::row_wise)
        , nz_(nz)
        , space_(0)
        , comm_(0)         
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

        clear();
      }

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
            std::cout << "Only one process, returning! \n";
            return;
          }

          // exchange data 
          DiscreteFunctionType tmp("ImprovedBCRSMatrix::communicate_tmp",*space_,arg);
          comm_->exchange( tmp );
        }
      }

      //! mult y = Ax 
      void mult(const BlockVectorType& x, BlockVectorType& y) const 
      {
        // exchange data 
        communicate( x );

        /*
        ConstRowIterator endi=this->end();
        for (ConstRowIterator i=this->begin(); i!=endi; ++i)
        {
          ConstColIterator endj = (*i).end();
          for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          {
            FMatrixHelp::multAssign((*j),x[j.index()],y[i.index()]);
          }
        }
        */

        // multiply 
        y = 0;
        this->umv(x,y);
      }
  };

  template <class RowSpaceType, class ColumnSpaceType> 
  class ISTLMatrixObject  
  {
    typedef typename RowSpaceType::GridType::template Codim<0>::Entity EntityType;
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

    const int numIterations_; 
    const double relaxFactor_; 
      
    const bool preconditioning_;

    //! setup matrix handler 
    ISTLMatrixObject(const RowSpaceType & rowSpace,
                     const ColumnSpaceType & colSpace,
                     bool preconditioning)
      : rowSpace_(rowSpace)
      , colSpace_(colSpace)
      , factor_((dim * 2) + 1)
      , size_(-1)
      , matrix_(0)
      , preconder_(0)
      , comm_(rowSpace_)
      , numIterations_(5)
      , relaxFactor_(1.1)
      , preconditioning_(preconditioning)
    {
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
        matrix_ = new MatrixType(rowSpace_, comm_, size_, colSpace_.indexSet().size(0),factor_);
        matrix().setup(rowSpace_,colSpace_);
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
