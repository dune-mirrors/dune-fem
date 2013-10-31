#ifndef DUNE_FEM_SPMATRIX_HH
#define DUNE_FEM_SPMATRIX_HH
 
//- system includes 
#include <vector>
#include <set>
#include <algorithm>

//- local includes 
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/operator/common/localmatrix.hh> 
#include <dune/fem/operator/common/localmatrixwrapper.hh> 
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#ifdef ENABLE_UMFPACK 
#include <umfpack.h>
// #include <ldl.h>
#endif

namespace Dune
{

  namespace Fem 
  {
  
    //! anonymous namespace, such that variable is only known within 
    //! SparseRowMatrix and other classes located here
    namespace
    {    
      //! If you have problems with this sparsematrix class, this might be due to
      //! inconsistencies produced in some methods.  
      //! In this case, you should turn on the consistencycheck of all non-const 
      //! methods by setting the following variable to 1 / 0 and by this locate the 
      //! buggy member method. Default is 0 = Check Off
      const int checkNonConstMethods = 0;
    }
    
    //*****************************************************************
    //
    //  --SparseRowMatrix
    //  
    //! Compressed row sparse matrix, where only the nonzeros of a row are
    //! keeped 
    //! (except if you "set" a single element explicitly 
    //! with the value 0, which is not forbidden and an element entry is 
    //! created)
    //!
    //*****************************************************************
      
    template <class T>
    class SparseRowMatrix
    {
    public: 
      typedef T Ttype;  //! remember the value type

      typedef SparseRowMatrix<T> ThisType;

      enum { defaultCol = -1 };
      enum { firstCol = defaultCol + 1 };

    protected:
      T* values_ ;      //! data values (nz_ * dim_[0] elements)
      int* col_;        //! row_ptr (dim_[0]+1 elements)
      int* nonZeros_;   //! row_ptr (dim_[0]+1 elements)
      int dim_[2];      //! dim_[0] x dim_[1] Matrix
      int nz_;          //! number of nonzeros per row

      int memSize_;
      bool sorted_;

      // temporary mem for resort 
      std::vector<int> newIndices_;
      std::vector<T> newValues_;

      std::set< int > clearedRows_;

      // omega for ssor preconditioning
      const double omega_;
      
      //! Copy Constructor prohibited 
      SparseRowMatrix(const SparseRowMatrix<T> &S);
    public:
      //! makes Matrix of zero length, omega is 1.1 by default 
      SparseRowMatrix(double omega = 1.1); 

      //! make matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row 
      //! and intialize all values with 'val' and set omega_ to omega
      SparseRowMatrix(int rows, int cols, int nz, const T& val = 0, 
                      double omega = 1.1);
     
      //! reserve memory for given rows, and number of non zeros, 
      //! set all entries to value dummy.... What is the use of this value?
      //! only initializing with 0 makes sense, so this argument is renamed
      //! 'dummy', by the way, nothing is happening with this argument, so 
      //! might be completely removed.
      void reserve(int rows, int cols, int nz, const T& dummy);

      //! resize keeping old values if possible, assuming rows == cols  
      void resize ( int newSize );

      //! resize keeping old values if possible   
      void resize ( int newRow, int newCol , int newNz = -1 );

      //! free memory for values_ and col_
      ~SparseRowMatrix();
        
    /*******************************/
    /*  Access and info functions  */
    /*******************************/
      //! length of array used to store matrix 
      int numberOfValues () const { return dim_[0]*nz_; }
      //! matrix value taken from real array 
      T&  val(int i) 
      { 
        assert( i >= 0 );
        assert( i < numberOfValues () );
        return values_[i]; 
      }

      //! return reference to value on given entry 
      const T&  val(int i) const { 
        assert( i >= 0 );
        assert( i < numberOfValues () );
        return values_[i]; 
      }

      //! return value and clear matrix entry 
      T popValue (int i) 
      { 
        if (checkNonConstMethods) assert(checkConsistency());
        T v = val(i); 
        values_[i] = 0;
        col_[i] = -1;
        if (checkNonConstMethods) assert(checkConsistency());
        return v;
      }

      //! return real column number for (row,localCol) 
      int realCol (int row, int fakeCol) const 
      {
        assert(fakeCol<dim_[1]);
        assert( row < dim_[0] );
        int pos = row*nz_ + fakeCol;
        return col_[pos];
      }

      //! return pair< value, column >, used by BlockMatrix 
      std::pair < T , int > realValue(int row, int fakeCol) 
      {
        assert( row < dim_[0] );
        assert( fakeCol < nz_ );
        int pos = row*nz_ + fakeCol;
        return realValue(pos);
      }
      
      //! return pair< value, column >, used by BlockMatrix 
      std::pair < const T , int > realValue(int row, int fakeCol) const 
      {
        assert( row < dim_[0] );
        assert( fakeCol < nz_ );
        int pos = row*nz_ + fakeCol;
        return realValue(pos);
      }
      
      //! return pair< value, column >, used by BlockMatrix 
      std::pair < T , int > realValue(int index) 
      {
        return std::pair< T, int > (values_[index], col_[index]); 
      }

      //! return pair< value, column >, used by BlockMatrix 
      std::pair < const T , int > realValue(int index) const 
      {
        return std::pair< const T, int > (values_[index], col_[index]); 
      }

      //! returns local col index for given global (row,col) 
      int colIndex(int row, int col);
      
      //! returns true if entry (row,col) exists in matrix 
      bool find (int row, int col) const;

      //! return number of rows = 0, cols = 1
      int dim(int i) const {return dim_[i];}
      //! return number of rows = 0, cols = 1
      int size(int i) const {return dim_[i];}
      
      //! return number of rows  
      int rows() const {return dim_[0];}

      //! return number of columns 
      int cols() const {return dim_[1];}

      //! return max number of non zeros 
      int numNonZeros() const {return nz_;}
      
      //! return number of non zeros in row 
      int numNonZeros(int i) const 
      { 
        assert( nonZeros_ );
        return nonZeros_[i]; 
      }

      //! return value of entry (row,col)
      T operator() ( const int row, const int col ) const;
      T operator() ( const unsigned int row, const unsigned int col ) const;
      T operator() ( const long unsigned int row, const long unsigned int col ) const
      {
        return this->operator()((unsigned int)(row), (unsigned int)(col) );
      }

      //! set entry to value
      //! note, that every entry is performed into the matrix!
      //! also setting of value 0 will result in an entry. So these
      //! calls should be ommited on a higher level 
      void set(int row, int col, T val);
      
      //! set all entries in row to zero 
      void clearRow (int row);

      //! ser all entris in column col to zero
      void clearCol ( int col );

      //! set all entries in row to zero 
      void scaleRow (int row, const T& val );

      //! set all matrix entries to zero, no other value makes sense for 
      //! sparse matrix 
      void clear();
      
      //! add value to row,col entry 
      void add(int row, int col, T val);
          
      //! muliply with scalar value 
      void multScalar(int row, int col, T val);
      
      //! make unitRow(row) and unitCol(col)
      void kroneckerKill(int row, int col);

      //! same as apply A * x = ret 
      template <class VECtype> 
      void mult(const VECtype *x, VECtype * ret) const;

      //! same as apply A * x = ret, used by OEM-Solvers 
      template <class VECtype> 
      T multOEMRow (const VECtype *x , const int row ) const;

      //! same as apply A * x = ret, used by OEM-Solvers 
      template <class VECtype> 
      void multOEM(const VECtype *x, VECtype * ret) const;

      //! calculates ret += A * x 
      template <class VECtype> 
      void multOEMAdd(const VECtype *x, VECtype * ret) const;

      //! same as apply A^T * x = ret, used by OEM-Solvers 
      template <class VECtype> 
      void multOEM_t(const VECtype *x, VECtype * ret) const;

      //! A(f) = ret, same as mult 
      template <class DiscFType, class DiscFuncType>
      void apply(const DiscFType &f, DiscFuncType &ret) const;
      
      //! A^T(f) = ret
      template <class ArgDFType, class DestDFType>
      void apply_t(const ArgDFType& f, DestDFType &ret) const;
      
      //! A(f) = ret 
      template <class DiscFuncType> 
      void operator () (const DiscFuncType &f, DiscFuncType &ret) const 
      {
        apply(f,ret); 
      }
      
      //! return diagonal of (this * A * B)
      template <class DiscFuncType>
      void getDiag(const ThisType&A, const ThisType &B, DiscFuncType &rhs) const;
      
      //! return diagonal of (this * A)
      template <class DiscFuncType>
      void getDiag(const ThisType &A, DiscFuncType &rhs) const;
      
      //! return diagonal entries of this matrix 
      template <class DiscFuncType>
      void getDiag(DiscFuncType &rhs) const;
      
      //! add diagonal to given DiscreteFunction
      template <class DiscFuncType>
      void addDiag(DiscFuncType &rhs) const;
      
      //! print matrix 
      void print (std::ostream& s) const;

      //! print values 
      void printReal (std::ostream& s) const;
         
      //! print columns
      void printColumns(std::ostream& s) const;
      
      //! print row-wise stored number of nonzeros 
      //! No counting is performed, but the member variable nonZeros_[]
      //! is reported. So here inconsistencies can occur to the true 
      //! nonzero entries in the matrix.
      void printNonZeros(std::ostream& s) const;
      
      //! check consistency, i.e. whether number of stored nonzeros
      //! corresponds to the counters in nonZeros[] and check, whether all
      //! columns within this range have reasonable values
      //! true == consistent
      //! false == non consistent
      //! an assert(checkConsistency()) can be called at entry and exit of
      //! non-const sparsematrix operations for ensuring maintaining of
      //! consistence. This can be made conditional by the member variable 
      //! checkNonConstMethods
      bool checkConsistency() const;

      //! make row a row with 1 on diagonal and all other entries 0 
      void unitRow(int row);

      //! make column a column with 1 on diagonal and all other entries 0 
      void unitCol(int col);

      //! check symetry 
      void checkSym ();

      // res = this * B 
      void multiply(const ThisType & B, ThisType & res) const;

      //! multiply this matrix with scalar 
      void scale(const T& factor);
      
      //! add other matrix to this matrix 
      void add(const ThisType & B ); 

      //! resort to have ascending column numbering 
      void resort();

      //! resort row to have ascending column numbering 
      void resortRow(const int row);

      //! SSOR preconditioning 
      void ssorPrecondition (const T*, T*) const;

      //! returns true if preconditioing is called before matrix multiply 
      bool rightPrecondition() const { return true; }
      
      //! apply preconditioning, calls ssorPreconditioning at the moment 
      void precondition (const T*u , T*x) const {  ssorPrecondition(u,x); }

      //! solve A x = b using the UMFPACK direct solver 
      void solveUMF(const T* b, T* x);
      void solveUMFNonSymmetric(const T* b, T* x);

    private:
      void setupUMF(int n, int nAll, int* Ap, int* Ai, T* Ax, int &ANZ, int &LNZ);
      //! delete memory 
      void removeObj();
    };


    // SparseRowMatrixObject
    // ---------------------
    template <class DomainSpace, class RangeSpace, 
              class Matrix = SparseRowMatrix< typename DomainSpace :: RangeFieldType > >
    class SparseRowMatrixObject;

    template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
    struct SparseRowMatrixTraits
    {
      typedef RowSpaceImp RowSpaceType;
      typedef ColSpaceImp ColumnSpaceType;
      typedef SparseRowMatrixTraits<RowSpaceType,ColumnSpaceType> ThisType;
      typedef SparseRowMatrixObject<RowSpaceType,ColumnSpaceType> MatrixObjectType;

      typedef RowSpaceImp  DomainSpaceType;
      typedef ColSpaceImp  RangeSpaceType;

    };

    template< class DomainSpace, class RangeSpace, class Matrix >
    class SparseRowMatrixObject
    : public Fem :: OEMMatrix
    {
    public:
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;

      /******************************************************************* 
      *   Rows belong to the DomainSpace and Columns to the RangeSpace   *
      *******************************************************************/    
      typedef typename DomainSpaceType :: EntityType  DomainEntityType ;
      typedef typename RangeSpaceType  :: EntityType  RangeEntityType ;
      typedef typename DomainSpaceType :: EntityType  RowEntityType ;
      typedef typename RangeSpaceType  :: EntityType  ColumnEntityType ;

      typedef typename DomainSpaceType :: BlockMapperType DomainBlockMapperType ;
      typedef NonBlockMapper< DomainBlockMapperType, 
                              DomainSpaceType :: localBlockSize > DomainMapperType;
      typedef typename RangeSpaceType :: BlockMapperType RangeBlockMapperType ;
      typedef NonBlockMapper< RangeBlockMapperType, 
                              RangeSpaceType :: localBlockSize > RangeMapperType;

    private:  
      typedef Matrix  MatrixType ;
      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, MatrixType > ThisType;

    protected:
      typedef typename DomainSpaceType :: GridType GridType;


      template< class MatrixObject >
      struct LocalMatrixTraits;

      template< class MatrixObject >
      class LocalMatrix;
      
    public:  
      typedef MatrixType PreconditionMatrixType;

    public:
      //! type of local matrix 
      typedef LocalMatrix<ThisType> ObjectType;
      typedef ThisType LocalMatrixFactoryType;
      typedef Fem :: ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
      //! type of local matrix 
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

      typedef ColumnObject< ThisType > LocalColumnObjectType;

    protected:
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;

      DomainMapperType domainMapper_ ;
      RangeMapperType  rangeMapper_ ;
      
      int sequence_;

      mutable MatrixType matrix_;
      bool preconditioning_;

      mutable LocalMatrixStackType localMatrixStack_;

    public:
      //! setup matrix handler 
      inline SparseRowMatrixObject( const DomainSpaceType &domainSpace,
                                    const RangeSpaceType &rangeSpace,
                                    const std::string &paramfile = "" )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        domainMapper_( domainSpace_.blockMapper() ),
        rangeMapper_( rangeSpace_.blockMapper() ),
        sequence_( -1 ),
        matrix_(),
        preconditioning_( false ),
        localMatrixStack_( *this )
      {
        int precon = 0;
        if( paramfile != "" )
        {
          readParameter( paramfile, "Preconditioning", precon );
        }
        else 
        {
          precon = Parameter :: getValue("Preconditioning", precon );
        }
        preconditioning_ = (precon > 0) ? true : false;
      }

      //! return domain space (i.e. space that builds the rows)
      const DomainSpaceType& domainSpace() const { return domainSpace_; }

      //! return range space (i.e. space that builds the columns)
      const RangeSpaceType& rangeSpace() const { return rangeSpace_; }

      //! return reference to stability matrix
      inline MatrixType &matrix () const
      {
        return matrix_;
      }

      //! interface method from LocalMatrixFactory
      inline ObjectType *newObject () const
      {
        return new ObjectType( *this, domainSpace_, rangeSpace_, domainMapper_, rangeMapper_ );
      }

      //! return local matrix 
      inline LocalMatrixType localMatrix( const DomainEntityType &domainEntity,
                                          const RangeEntityType &rangeEntity ) const
      {
        /******************************************************************* 
        *   Rows belong to the DomainSpace and Columns to the RangeSpace   *
        *******************************************************************/
        return LocalMatrixType( localMatrixStack_, domainEntity, rangeEntity );
      }

      LocalColumnObjectType localColumn( const DomainEntityType &domainEntity ) const
      {
        return LocalColumnObjectType ( *this, domainEntity );
      }

      //! resize all matrices and clear them 
      inline void clear ()
      {
        matrix_.clear();
      }

      //! return true if precoditioning matrix is provided 
      bool hasPreconditionMatrix () const
      {
        return preconditioning_;
      }

      //! return reference to preconditioner 
      const PreconditionMatrixType &preconditionMatrix () const 
      { 
        return matrix_;
      }

      //! reserve memory for assemble based on the provided stencil 
      template <class Stencil>
      inline void reserve(const Stencil &stencil, bool verbose = false )
      {
        if( sequence_ != domainSpace_.sequence() )
        {
#ifndef DUNE_FEM_DONT_CHECKITERATORS_OF_SPACE
          // if empty grid do nothing (can appear in parallel runs)
          if( (domainSpace_.begin() != domainSpace_.end())
              && (rangeSpace_.begin() != rangeSpace_.end()) )
#endif
          {
            
            if( verbose )
            {
              const int rowMaxNumbers = rangeMapper_.maxNumDofs();
              const int colMaxNumbers = domainMapper_.maxNumDofs();

              std::cout << "Reserve Matrix with (" << rangeSpace_.size() << "," << domainSpace_.size()<< ")\n";
              std::cout << "Max number of base functions = (" << rowMaxNumbers << "," << colMaxNumbers << ")\n";
            }

            // upper estimate for number of non-zeros 
            const static size_t domainLocalBlockSize = DomainSpaceType::localBlockSize;
            const int nonZeros = std::max( (int)(stencil.maxNonZerosEstimate()*domainLocalBlockSize), 
                                           matrix_.numNonZeros() );
            matrix_.reserve( rangeSpace_.size(), domainSpace_.size(), nonZeros, 0.0 );
          }
          sequence_ = domainSpace_.sequence();
        }
      }

      //! solve A dest = arg using the UMFPACK direct solver 
      template< class DomainFunction, class RangeFunction >
      void solveUMF( const DomainFunction &arg, RangeFunction &dest ) const
      {
        DUNE_THROW(NotImplemented,"solveUMF only implemented for AdaptiveDiscreteFunctions");
      }
      template< class DomainFunction, class RangeFunction >
      void solveUMFNonSymmetric( const DomainFunction &arg, RangeFunction &dest ) const
      {
        DUNE_THROW(NotImplemented,"solveUMF only implemented for AdaptiveDiscreteFunctions");
      }

      //! solve A dest = arg using the UMFPACK direct solver 
      void solveUMF ( const AdaptiveDiscreteFunction< DomainSpaceType > &arg, 
                      AdaptiveDiscreteFunction< RangeSpaceType> &dest ) const
      {
        matrix_.solveUMF( arg.leakPointer(), dest.leakPointer() );
      }
      //! solve A dest = arg using the UMFPACK direct solver 
      void solveUMFNonSymmetric ( const AdaptiveDiscreteFunction< DomainSpaceType > &arg, 
                                  AdaptiveDiscreteFunction< RangeSpaceType> &dest ) const
      {
        matrix_.solveUMFNonSymmetric( arg.leakPointer(), dest.leakPointer() );
      }

      //! apply matrix to discrete function
      template< class DomainFunction, class RangeFunction >
      void apply ( const DomainFunction &arg, RangeFunction &dest ) const
      {
        // do matrix vector multiplication 
        matrix_.apply( arg, dest );

        // communicate data 
        dest.communicate();
      }

      //! apply matrix to discrete function
      void apply ( const AdaptiveDiscreteFunction< DomainSpaceType > &arg, 
                   AdaptiveDiscreteFunction< RangeSpaceType> &dest ) const
      {
        // do matrix vector multiplication 
        matrix_.multOEM( arg.leakPointer(), dest.leakPointer() );

        // communicate data 
        dest.communicate();
      }

      //! apply transposed matrix to discrete function
      template< class DomainFunction, class RangeFunction >
      void apply_t ( const RangeFunction &arg, DomainFunction &dest ) const
      {
        // do matrix vector multiplication 
        matrix_.apply_t( arg, dest );

        // communicate data 
        dest.communicate();
      }

      //! apply transposed matrix to discrete function
      void apply_t ( const AdaptiveDiscreteFunction< RangeSpaceType > &arg, 
                     AdaptiveDiscreteFunction< DomainSpaceType> &dest ) const
      {
        // do matrix vector multiplication 
        matrix_.multOEM_t( arg.leakPointer(), dest.leakPointer() );

        // communicate data 
        dest.communicate();
      }


      //! mult method of matrix object used by oem solver
      double ddotOEM( const double *v, const double *w ) const
      {
        typedef AdaptiveDiscreteFunction< DomainSpaceType > DomainFunctionType;
        DomainFunctionType V( "ddot V", domainSpace_, v );
        DomainFunctionType W( "ddot W", domainSpace_, w );
        return V.scalarProductDofs( W );
      }

      //! mult method of matrix object used by oem solver
      void multOEM( const double *arg, double *dest ) const
      {
        typedef AdaptiveDiscreteFunction< DomainSpaceType > DomainFunctionType;
        typedef AdaptiveDiscreteFunction< RangeSpaceType > RangeFunctionType;

        DomainFunctionType farg( "multOEM arg", domainSpace_, arg );
        RangeFunctionType fdest( "multOEM dest", rangeSpace_, dest );
        apply( farg, fdest );
      }

      //! resort row numbering in matrix to have ascending numbering 
      void resort() 
      {
        matrix_.resort();
      }

      void createPreconditionMatrix()
      { 
      }

      //! extract diagonal entries from matrix into discrete function 
      template < class DiscreteFunctionType > 
      void extractDiagonal( DiscreteFunctionType& diag ) const
      {
        // this only works for matrices with same number of rows,cols
        assert( matrix_.rows() == matrix_.cols() );
        typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType ;
        const DofIteratorType dofEnd = diag.dend();
        unsigned int row = 0;
        for( DofIteratorType dofIt = diag.dbegin(); 
             dofIt != dofEnd; ++ dofIt, ++row ) 
        {
          assert( row < ( unsigned int )matrix_.rows() );
          (*dofIt) = matrix_( row, row );
        }
      }

      /** \brief delete all row belonging to a hanging node and rebuild them */
      template <class HangingNodesType> 
      void changeHangingNodes(const HangingNodesType& hangingNodes) 
      {
        {
          typedef typename HangingNodesType :: IteratorType IteratorType;
          const IteratorType end = hangingNodes.end();
          for( IteratorType it = hangingNodes.begin(); it != end; ++it)
          {
            insertHangingRow( hangingNodes, (*it).first , (*it).second );
          }
        }

        /*
        {
          typedef typename HangingNodesType :: SlaveNodesType SlaveNodesType; 
          typedef typename SlaveNodesType :: const_iterator iterator;
          iterator end =  hangingNodes.slaveNodes().end();
          for( iterator it =  hangingNodes.slaveNodes().begin(); 
               it != end; ++it )
          {
            matrix().unitRow( *it );
            matrix().set( *it, *it, 0.0 );
          }
        }
        */
      }
  protected:
      /** \brief insert row to be a row for a hanging node */
      template <class HangingNodesType, class ColumnVectorType>
      void insertHangingRow( const HangingNodesType& hangingNodes,
                             const int row, const ColumnVectorType& colVec)
      {
        const size_t cols = colVec.size();

        // distribute row to associated rows 
        const int nonZeros = matrix().numNonZeros( row );
        for( int c = 0; c < nonZeros; ++c)
        {
          std::pair< double, int > val =  matrix().realValue( row, c );
          for( size_t j = 0; j < cols; ++ j)
          {
            const double value = colVec[j].second * val.first;
            assert( ! hangingNodes.isHangingNode( colVec[j].first ) );
            matrix().add( colVec[j].first , val.second, value );
          }
        }

        // replace hanging row 
        matrix().unitRow( row );
        for( size_t j = 0; j < cols; ++ j)
        {
          matrix().add( row, colVec[j].first , -colVec[j].second );
        }
      }
    };



    template< class DomainSpace, class RangeSpace, class Matrix >
    template< class MatrixObject >
    struct SparseRowMatrixObject< DomainSpace, RangeSpace, Matrix >::LocalMatrixTraits
    {
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;

      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, Matrix >
        SparseRowMatrixObjectType;
      
      typedef typename SparseRowMatrixObjectType :: template LocalMatrix< MatrixObject > LocalMatrixType;

      typedef typename RangeSpaceType :: RangeFieldType RangeFieldType;
      typedef RangeFieldType LittleBlockType;

      typedef typename SparseRowMatrixObjectType::DomainMapperType  DomainMapperType;
      typedef typename SparseRowMatrixObjectType::RangeMapperType   RangeMapperType;
    };



    //! LocalMatrix 
    template< class DomainSpace, class RangeSpace, class Matrix >
    template< class MatrixObject >
    class SparseRowMatrixObject< DomainSpace, RangeSpace, Matrix > :: LocalMatrix
    : public LocalMatrixDefault< LocalMatrixTraits< MatrixObject > >
    {
    public:
      //! type of matrix object 
      typedef MatrixObject MatrixObjectType;

      //! type of the traits
      typedef LocalMatrixTraits< MatrixObjectType > Traits;

    private:
      typedef LocalMatrixDefault< Traits > BaseType;

    public:
      //! type of matrix 
      typedef typename MatrixObjectType :: MatrixType MatrixType;

      //! type of entries of little blocks 
      typedef typename Traits :: RangeFieldType RangeFieldType;

      //! type of the DoFs
      typedef RangeFieldType DofType;

      //! type of little blocks 
      typedef typename Traits :: LittleBlockType LittleBlockType;

      //! type of nonblocked domain mapper   
      typedef typename Traits :: DomainMapperType DomainMapperType;
      //! type of nonblocked domain mapper   
      typedef typename Traits :: RangeMapperType RangeMapperType;

    protected:
      MatrixType &matrix_; 
      const DomainMapperType& domainMapper_;
      const RangeMapperType&  rangeMapper_;

      typedef std :: vector< typename RangeMapperType :: SizeType > RowIndicesType ;
      //! global index in the DomainSpace
      RowIndicesType rowIndices_;

      typedef std :: vector< typename DomainMapperType :: SizeType > ColumnIndicesType ;
      //! global index in the RangeSpace
      ColumnIndicesType columnIndices_;

      using BaseType :: domainSpace_;
      using BaseType :: rangeSpace_;
      
    public:  
      //! constructor taking entity and spaces for using mapToGlobal
      //! class RowSpaceType, class ColSpaceType> 
      inline LocalMatrix( const MatrixObjectType &matrixObject,
                          const DomainSpaceType &domainSpace,
                          const RangeSpaceType &rangeSpace,
                          const DomainMapperType& domainMapper,
                          const RangeMapperType& rangeMapper )
      : BaseType( domainSpace, rangeSpace),
        matrix_( matrixObject.matrix() ),
        domainMapper_( domainMapper ),
        rangeMapper_( rangeMapper )
      {
      }
      
    private: 
      // prohibit copying 
      LocalMatrix( const LocalMatrix & );

    public:
      void init( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        /******************************************************************* 
        *   Rows belong to the DomainSpace and Columns to the RangeSpace   *
        *******************************************************************/

        // initialize base functions sets 
        BaseType::init( domainEntity, rangeEntity );
          
        // rows are determined by the range space
        rowIndices_.resize( rangeMapper_.numDofs( rangeEntity ) );
        rangeMapper_.mapEach( rangeEntity, AssignFunctor< RowIndicesType >( rowIndices_ ) );

        // columns are determind by the domain space
        columnIndices_.resize( domainMapper_.numDofs( domainEntity ) );
        domainMapper_.mapEach( domainEntity, AssignFunctor< ColumnIndicesType >( columnIndices_ ) );
      }

      //! return number of rows 
      int rows () const
      {
        return rowIndices_.size();
      }

      //! return number of columns 
      int columns () const
      {
        return columnIndices_.size();
      }

      //! add value to matrix entry
      void add( int localRow, int localCol, const DofType value )
      {
        assert( value == value );
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        matrix_.add( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! get matrix entry 
      DofType get( int localRow, int localCol ) const
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        return matrix_( rowIndices_[ localRow ], columnIndices_[ localCol ] );
      }

      //! set matrix entry to value 
      void set( int localRow, int localCol, const DofType value )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        matrix_.set( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! set matrix row to zero except diagonla entry 
      void unitRow( const int localRow )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.unitRow( rowIndices_[ localRow ] );
      }

      //! set matrix row to zero
      void clearRow( const int localRow )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.clearRow( rowIndices_[localRow]);
      }

      //! set matrix column to zero
      void clearCol ( const int localCol )
      {
        assert( (localCol >= 0) && (localCol < columns()) );
        matrix_.clearCol( columnIndices_[localCol] );
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        const int row = rows();
        for( int i = 0; i < row; ++i )
          matrix_.clearRow( rowIndices_[ i ] );
      }

      //! scale local matrix with a certain value 
      void scale ( const DofType& value ) 
      {
        const int row = rows();
        for( int i = 0; i < row; ++i )
          matrix_.scaleRow( rowIndices_[ i ] , value );
      }

      //! resort all global rows of matrix to have ascending numbering 
      void resort ()
      {
        const int row = rows();
        for( int i = 0; i < row; ++i )
          matrix_.resortRow( rowIndices_[ i ] );
      }
    };

  } // namespace Fem

} // namespace Dune

#include "spmatrix.cc"

#endif // #ifndef DUNE_FEM_SPMATRIX_HH
