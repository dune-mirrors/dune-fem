#ifndef DUNE_FEOPERATOR_HH
#define DUNE_FEOPERATOR_HH

//- Dune includes 
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/referenceelements.hh>

//- local includes 
#include "common/operator.hh"
#include "common/localoperator.hh"
#include "feop/spmatrix.hh"


namespace Dune {

/** @defgroup FEOpInterface FEOpInterface 
    @ingroup DiscreteOperator
 
   
  @{
 */

/*======================================================================*/
/*!
 *  \class FEOpInterface
 *  \brief FEopInterface is the interface for the definition of a finite 
 *         element operator.
 */
/*======================================================================*/


template <class DiscFunctionType, class FEOpImp>
class FEOpInterface 
: public Operator<typename DiscFunctionType::DomainFieldType, 
    typename DiscFunctionType::RangeFieldType,DiscFunctionType,DiscFunctionType> 
{
public:

/*======================================================================*/
/*! 
 *   getLocalMatrix: Interface method that returns the local matrix of the 
 *                   finite element operator on an entity
 *
 *   This method has to be provided by derived classes.
 *
 *   \param the entity, the local matrix size and a reference to an 
 *          instance of the local matrix implementation
 */
/*======================================================================*/

  template <class  EntityType, class LocalMatrixImp>
  void getLocalMatrix( EntityType &entity, const int matSize, LocalMatrixImp& mat) const {
    return asImp().getLocalMatrix( entity, matSize, mat);
  }

protected:
  // Barton-Nackman
  FEOpImp &asImp() { return static_cast<FEOpImp&>( *this ); }

  const FEOpImp &asImp( ) const { return static_cast<const FEOpImp&>( *this ); }
};


/** @} end documentation group */

/*======================================================================*/
/*!
 *  \class FEOp 
 *  \brief The FEOp class provides one example of a class satisfying the 
 *         FEOpInterface. 
 *  
 *  The MatrixImp class must provide a storage type for the global 
 *  operator matrix and some arithmetics and basic functionality. In 
 *  particular a print() method and an apply() method are assumed to exist.
 *  Additionally a constructor with three arguments is required as 
 *  explained in newEmptyMatrix. An add(row,col,val) method is assumed 
 *  to exist.
 *
 *  different operating modes are possible. In case of ASSEMBLED, the 
 *  whole 
 *  global operator matrix is allocated and completely precomputed by 
 *  corresponding methods. In case of ON_THE_FLY, no complete matrix is
 *  allocated, but the matrix-vector multiplication is performed by 
 *  on-the-fly computation of the elementwise local matrices.
 */
/*======================================================================*/

template <class DiscFunctionType, class MatrixImp, class FEOpImp>
class FEOp : public FEOpInterface<DiscFunctionType,FEOpImp> ,
public LocalOperatorDefault <DiscFunctionType,DiscFunctionType, typename
DiscFunctionType::RangeFieldType , FEOpImp  >
{
  
public:
  //! Type of matrix storage class used for global operator matrix
  typedef MatrixImp MatrixType;  
  
  //! Operation mode: global allocation and multiplication or only 
  //! on-the-fly
  enum OpMode { ON_THE_FLY, ASSEMBLED };

  //! fix the size of the local matrices to some reasonable extent
  enum { maxnumOfBaseFct = 100 };
  
/*======================================================================*/
/*! 
 *   constructor: Initialization of FEOp
 *
 *   Based on an existing instance of a function space the FEOp is 
 *   initialized. 
 *   Operation mode must be selected as ASSEMBLED or ON_THE_FLY. 
 *   
 *   ???? The role of isleaf is unclear ?????
 * 
 *   \param an instance of the discrete function space, the operator mode 
 *          and a leaf-flag
 *
 *   \return the initialized FEOp
 */
/*======================================================================*/

  FEOp( const typename DiscFunctionType::FunctionSpaceType &fuspace, 
                         OpMode opMode = ASSEMBLED, bool leaf = true ) :
    functionSpace_( fuspace ),  matrix_ (0), matrix_assembled_( false ),
    arg_ ( NULL ), dest_ (NULL) , opMode_(opMode) , leaf_ (leaf) {};

  //! \todo Please doc me!
/*======================================================================*/
/*! 
 *   destructor: In case of allocation of global operator matrix
 *               it is deallocated.
 */
/*======================================================================*/

  ~FEOp( ) {
    if ( matrix_ ) delete matrix_;
  };

public:
/*======================================================================*/
/*! 
 *   print: print matrix to standard out 
 */
/*======================================================================*/

  void print () const 
  {
    if(!this->matrix_assembled_) this->assemble();
    this->matrix_->print(std::cout);
  }
 
/*======================================================================*/
/*! 
 *   myMatrix: return reference to Matrix for oem solvers. 
 *
 *   The assembled matrix is returned. That means, if global matrix is not 
 *   yet allocated, a new empty matrix is generated. If current global 
 *   matrix is not yet assembled, an assembly is initiated.
 *
 *   ??? What happens in case of ON_THE_FLY ??? 
 *
 *   \return reference to assembled global matrix
 */
/*======================================================================*/

  MatrixType & myMatrix() const
  {
    //assert(matrix_assembled_ == true);
    if ( !this->matrix_assembled_ )
    {
      if(!this->matrix_)
        this->matrix_ = this->newEmptyMatrix( );
      this->assemble(); 
    }
    return (*this->matrix_);
  }
 

/*======================================================================*/
/*! 
 *   hasPcMatrix: return false
 *
 *   ???  What is this method good for ???
 *
 *   \return false
 */
/*======================================================================*/

  bool hasPcMatrix() const { return false; }


/*======================================================================*/
/*! 
 *   pcMatrix: call myMatrix
 *              
 *   ???? The use and relevance of this method is unclear ????
 *
 *   \return the result of myMatrix
 */
/*======================================================================*/

  MatrixType & pcMatrix() const { return myMatrix(); }


/*======================================================================*/
/*! 
 *   initialize: delete and deallocate global matrix
 *
 *   ??? Why is this method called init, if it does the same as a 
 *       destructor should do ??? 
 *
 */
/*======================================================================*/
  
  //! methods for global application of the operator
  void initialize(){
    //std::cout << "Matrix reinitialized!" << std::endl ;
    
    matrix_assembled_ = false;
    if ( matrix_ ) delete(matrix_);
    matrix_ = 0;
  };


/*======================================================================*/
/*! 
 *   operator(): application operator
 *
 *   This method only makes sense in case of ASSEMBLED opmode, as a global
 *   apply() on the matrix is called, i.e. a matrix-vector-multiplication 
 *   with the vector arg is performed, the result stored in dest. This is 
 *   an implicit requirement on the MatrixImp class!
 *
 *   \param the argument for the matrix-vector multiplication and the 
 *          storage for the destination vector
 */
/*======================================================================*/

  virtual void operator()(const DiscFunctionType &arg, 
                          DiscFunctionType &dest) const 
  {
    assert( this->opMode_ == ASSEMBLED );
   
    if ( !matrix_assembled_ ) 
      {
        assemble(); 
      }

    matrix_->apply( arg, dest );
  };


public:
/*======================================================================*/
/*! 
 *   isLeaf: returns true if Leafiterator should be used, else false is 
 *   returned
 *
 *   ??? what is this good for? ???
 *
 *   \return the current isleaf_ value
 */
/*======================================================================*/

  bool isLeaf (){
    return leaf_;
  };

/*======================================================================*/
/*! 
 *  prepareGlobal: store argument and destination 
 *
 *  The storage for argument and destination is asserted and saved and
 *  the destination vector is cleared.
 *
 *   \param argument and destination storage 
 */
/*======================================================================*/
  void prepareGlobal(const DiscFunctionType &Arg, DiscFunctionType & Dest )  
  { 
    arg_  = &Arg.argument();
    dest_ = &Dest.destination();
    assert(arg_ != NULL); assert(dest_ != NULL);
    dest_->clear();
  };

  //! 

/*======================================================================*/
/*! 
 *   finalizeGlobal: detach the argument and destination storage
 * 
 *   The members of the argument and dest-pointer are set to NULL hereby
 *   detaching the user specified storage from any future operations of 
 *   FEOp. 
 */
/*======================================================================*/

  void finalizeGlobal()  
  { 
    arg_ = NULL; dest_ = NULL;
  };

  //! 
/*======================================================================*/
/*! 
 *   applyLocal: makes local multiply on the fly
 *
 *   the global operator matrix is not set up, but the local element 
 *   matrices are used for incremental matrix-vector multiplication.
 *   A grid-walkthrough and applyLocal results in a complete matrix-vector
 *   multiplication
 * 
 *   \param the entity
 */
/*======================================================================*/

  template <class EntityType>
  void applyLocal ( EntityType &en ) const 
  {
    const DiscFunctionType & arg  = (*arg_);
    DiscFunctionType & dest = (*dest_);

    typedef typename DiscFunctionType::FunctionSpace FunctionSpaceType;
    typedef typename FunctionSpaceType::GridType GridType; 
    
    typedef typename EntityType::IntersectionIterator NeighIt;
    
    typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    typedef typename FunctionSpaceType::RangeType RangeVecType;
    typedef typename FunctionSpaceType::JacobianRange JacobianRange;
    typedef typename FunctionSpaceType::DomainType DomainVecType;

    typedef typename DiscFunctionType::DofIteratorType DofIteratorType;
    typedef typename DiscFunctionType::ConstDofIteratorType ConstDofIteratorType;

    DofIteratorType dest_it = dest.dbegin();
    ConstDofIteratorType arg_it = arg.dbegin();
      
    const BaseFunctionSetType & baseSet = functionSpace_.getBaseFunctionSet( en );
    int numOfBaseFct = baseSet.numBaseFunctions();  

    assert( numOfBaseFct <= maxnumOfBaseFct );

    FieldMatrix<double, maxnumOfBaseFct, maxnumOfBaseFct> mat;
    
    getLocalMatrix( en, numOfBaseFct, mat);

    if(this->scalar_ == 1.)
    {
      for(int i=0; i<numOfBaseFct; i++) 
      {  
        int row = functionSpace_.mapToGlobal( en , i );
        for (int j=0; j<numOfBaseFct; j++ ) 
        {
          int col = functionSpace_.mapToGlobal( en , j );   

          // scalar comes from LocalOperatorDefault, if operator is scaled,
          // i.e. with timestepsize
          dest_it[ row ] += arg_it[ col ] * mat[i][j];
        }
      }
    }
    else 
    {
      for(int i=0; i<numOfBaseFct; i++) 
      {  
        int row = functionSpace_.mapToGlobal( en , i );
        for (int j=0; j<numOfBaseFct; j++ ) 
        {
          int col = functionSpace_.mapToGlobal( en , j );   

          // scalar comes from LocalOperatorDefault, if operator is scaled,
          // i.e. with timestepsize
          double val = (this->scalar_) * mat[i][j];

          dest_it[ row ] += arg_it[ col ] * val;
        }
      }
    }
  }; // end applyLocal



/*======================================================================*/
/*! 
 *   finalizeLocal: corrects the mapping in order to take into account 
 *                  dirichlet boundary conditions.
 * 
 *   The dofs are assumed to be correlated to the local vertex numbers of 
 *   the element, so only Lagrange basis are possible by this.
 *   The dirichlet-treatment is then simple copying of the argument DOF 
 *   to the destination DOF in case of Dirichlet boundary points. So the 
 *   assumption is, that the argument already contains correct 
 *   Dirichlet values.
 *
 *   \param the entity on which correction is to be performed
 */
/*======================================================================*/

 
  template <class EntityType>
  void finalizeLocal ( EntityType &en ) const 
  {
    // eliminate the Dirichlet rows and columns 
    typedef typename DiscFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
 
    const DiscFunctionType & arg  = (*arg_);
    DiscFunctionType & dest = (*dest_);

    const GridPartType & gridPart = arg.getFunctionsSpace().gridPart();

    typedef typename DiscFunctionType::DofIteratorType DofIteratorType;
    typedef typename DiscFunctionType::ConstDofIteratorType ConstDofIteratorType;

    DofIteratorType dest_it = dest.dbegin();
    ConstDofIteratorType arg_it = arg.dbegin();

    const GeometryType t = en.geometry().type();

    const IntersectionIteratorType endnit = gridPart.iend(en);
    for(IntersectionIteratorType nit = gridPart.ibegin(en); nit != endnit ; ++nit)
    {
      if(nit.boundary())
      {
        int face = nit.numberInSelf();
        enum { dim = EntityType :: dimension };
        typedef typename EntityType :: ctype coordType;
        
        const int faceCodim = 1;

        if(t.isSimplex())
        {
          if( nit.boundaryId() != 0 )
          {
            static ReferenceSimplex< coordType, dim > refElem;
            int novx = refElem.size( face, faceCodim , dim );
            assert( novx == dim );
            for(int j=0; j<novx ; j++)
            {
              // get all local numbers located on the face 
              int vx  = refElem.subEntity(face, faceCodim , j , dim );
              // get global dof numbers of this vertices 
              int col = functionSpace_.mapToGlobal( en, vx);
              // set solution on dirichlet bnd 
              dest_it[col] = arg_it[col];
            }
          }
        }
        if(t.isCube())
        {
          static ReferenceCube< coordType, dim > refElem;
          int novx = refElem.size( face, faceCodim , dim );
          for(int j=0; j<novx ; j++)
          {
            // get all local numbers located on the face 
            int vx  = refElem.subEntity(face, faceCodim , j , dim );
            // get global dof numbers of this vertices 
            int col = functionSpace_.mapToGlobal( en, vx );
            // set solution on dirichlet bnd 
            dest_it[col] = arg_it[col];
          }
        }
      }
    }
  }// end finalizeLocal

protected:
  //! the corresponding function_space 
  const typename DiscFunctionType::FunctionSpaceType & functionSpace_;
  
  //! pointer to the representing global matrix 
  mutable MatrixType *matrix_ ;
 
  //! flag indicating whether the global matrix is assembled 
  mutable bool matrix_assembled_;
 
  //! pointers to storage of argument and destination
  const DiscFunctionType * arg_;
  DiscFunctionType * dest_;

/*======================================================================*/
/*! 
 *   newEmptyMatrix: allocation of a new global matrix
 *
 *   The Matrixclass is required to have a constructor with the syntax
 *   MatrixType(nrows, ncols, nonzeros_per_row).   
 *
 *   \return a pointer to the newly allocated global matrix.
 */
/*======================================================================*/

  MatrixType* newEmptyMatrix( ) const 
  { typedef typename DiscFunctionType::FunctionSpaceType::GridType GridType; 
    enum { dim = GridType::dimension };
    return new MatrixType( this->functionSpace_.size ( ) , 
        this->functionSpace_.size ( ) , 
        15 * (dim-1));
  };

/*======================================================================*/
/*! 
 *   assemble: perform grid-walkthrough and assemble global matrix
 * 
 *   If the matrix storage is 
 *   not allocated, new storage is allocated by newEmptyMatrix.
 *   the begin and end iterators are determined and the assembling
 *   of the global matrix initiated by call of assembleOnGrid and 
 *   bndCorrectOnGrid. The assemled flag is set. 
 */
/*======================================================================*/

  void assemble ( ) const 
  {
    if(!this->matrix_) matrix_ = this->newEmptyMatrix( );
    
    typedef typename DiscFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::GridType GridType; 
    typedef typename FunctionSpaceType::IteratorType IteratorType;
   
    {  
      FieldMatrix<double, maxnumOfBaseFct, maxnumOfBaseFct> mat;

      IteratorType it    = functionSpace_.begin(); 
      IteratorType endit = functionSpace_.end(); 

      assembleOnGrid(it, endit, mat);
    }

    {
      IteratorType it    = functionSpace_.begin(); 
      IteratorType endit = functionSpace_.end(); 
      bndCorrectOnGrid(it, endit);
    }

    matrix_assembled_ = true;
  };
  
/*======================================================================*/
/*! 
 *   assembleOnGrid: perform grid walkthrough and assemble matrix
 *
 *   For each element, the local element matrix is determined into the 
 *   given local matrix storage and distributed into the global matrix.
 *   Distribution is performed by an add(row,col,val) method on the 
 *   global matrix class.
 *
 *   ??? why shoudl this one not be private? method assemble() is 
 *   sufficient for public call. ???
 *
 *   \param start and end iterator and storage for a local matrix 
 */
/*======================================================================*/

  template <class GridIteratorType, class LocalMatrixImp>
  void assembleOnGrid ( GridIteratorType &it, GridIteratorType &endit, 
                        LocalMatrixImp &mat) const
        {
          typedef typename DiscFunctionType::FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    // run through grid and add up local contributions
    for( ; it != endit; ++it )
    {
      const BaseFunctionSetType & baseSet = functionSpace_.getBaseFunctionSet( *it );
      const int numOfBaseFct = baseSet.numBaseFunctions();  
      
      // setup matrix 
      getLocalMatrix( *it, numOfBaseFct, mat);

      for(int i=0; i<numOfBaseFct; i++) 
      { 
        int row = functionSpace_.mapToGlobal( *it , i );
        for (int j=0; j<numOfBaseFct; j++ ) 
        {
          int col = functionSpace_.mapToGlobal( *it , j );    
          matrix_->add( row , col , mat[i][j]);
        }
      }
    }
  }

/*======================================================================*/
/*! 
 *   bndCorrectOnGrid: treatment of Dirichlet-DOFS
 *
 *   delete rows and columns for dirichlet DOFS, setting diagonal 
 *   element to 1. Lagrange Basis is implicitly assumed.
 *
 *   \param start and end iterator
 */
/*======================================================================*/

  
  //! \todo Please doc me!
  template <class GridIteratorType>
  void bndCorrectOnGrid ( GridIteratorType &it, const GridIteratorType &endit) const
  {
    // eliminate the Dirichlet rows and columns 
    typedef typename DiscFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::GridType GridType; 
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename FunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
 
    const GridPartType & gridPart = functionSpace_.gridPart();

    for( ; it != endit; ++it ) 
    {
      const EntityType & en = *it; 

      const GeometryType t = en.geometry().type();
      const IntersectionIteratorType endnit = gridPart.iend(en);
      for(IntersectionIteratorType nit = gridPart.ibegin(en); nit != endnit ; ++nit)
      {
        if(nit.boundary())
        {
          const int faceCodim = 1;
          int face = nit.numberInSelf();
          enum { dim = EntityType :: dimension };
          typedef typename EntityType :: ctype coordType;

          if( t.isSimplex() )
          {
            if( nit.boundaryId() != 0 )
            {
              static ReferenceSimplex< coordType, dim > refElem;
              int novx = refElem.size( face, faceCodim , dim );
              assert( novx == dim );
              for(int j=0; j<novx ; j++)
              {
                // get all local numbers located on the face 
                int vx  = refElem.subEntity(face, faceCodim , j , dim );
                // get global dof numbers of this vertices 
                int col = functionSpace_.mapToGlobal( en, vx);
                // set solution on dirichlet bnd 

                // unitRow unitCol for boundary
                matrix_->kroneckerKill(col,col);
              }
            }
          }
          if( t.isCube() )
          {
            static ReferenceCube< coordType, dim > refElem;
            int novx = refElem.size( face, faceCodim , dim );
            for(int j=0; j<novx ; j++)
            {
              // get all local numbers located on the face 
              int vx  = refElem.subEntity(face, faceCodim , j , dim );
              // get global dof numbers of this vertices 
              int col = functionSpace_.mapToGlobal( en, vx);

              // unitRow unitCol for boundary
              matrix_->kroneckerKill(col,col);
            }
          }
        }
      }
    }
  };

private:
  //! operator mode 
  OpMode opMode_;

  //! true if LeafIterator is used, deprecated because the Iterator now
  //! comes from the space 
  bool leaf_;
};


} // end namespace


#endif
