#ifndef DUNE_FEOP_HH
#define DUNE_FEOP_HH

//- Dune includes 
#include <dune/common/fmatrix.hh>

#include <dune/fem/operator/common/operator.hh>
// #include <dune/fem/operator/common/localoperator.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>
// #include<dune/fem/io/file/ioutils.hh>

// include temporary old and new sparsematrix class
#include <dune/fem/operator/matrix/spmatrix.hh>

namespace Dune
{

  /** \class FEOp 
   *  \ingroup EllipticOperator
   *  \brief The FEOp class provides an example of a Finite Element Operator 
   *     
   *  The class is an operator used for general elliptic problems, specialized 
   *  to Lagrange-bases, as explicit boundary setting and kronecker-kills
   *  in the matrix are performed.
   *
   *  The derivation from the FEOpInterface seems superfluous, as the local 
   *  element matrix provider is now specified by a template-parameter. 
   *
   *  The class is used for general elliptic problems + boundary treatment: 
   *  \f{displaymath}
   *  \begin{array}{rcll}
   *  -\nabla \cdot (a \nabla u + b u) + c u &=& f   & \quad \mbox{in} \enspace \Omega\\
   *                                       u &=& g_D & \quad \mbox{in} \enspace \Gamma_D\\
   *                    (a \nabla u + b u) n &=& g_N & \quad \mbox{in} \enspace \Gamma_N\\
   *         (a \nabla u + b u) n + \alpha u &=& g_R & \quad \mbox{in} \enspace \Gamma_R 
   *  \end{array}
   *  \f}
   *  where \f$a, b, c, g_D, g_N, g_R\f$ are space dependent, \f$\alpha\f$ a constant
   *  and the quantities denote
   *  \f$a( x )\f$:      stiffness
   *  \f$b( x )\f$:      velocity
   *  \f$c( x )\f$:      mass
   *  \f$f( x )\f$:      source
   *  \f$g_D\f$:         Dirichlet values
   *  \f$g_N\f$:         Neumann values
   *  \f$g_R\f$:         Robin values
   *  \f$\alpha( x )\f$: Robin alpha
   *
   * the following assumptions on the basis/model are made:
   *    - The discrete function space is a nodal basis, 
   *      there exists a set of x_i such that phi_j(x_i) = kronecker(i,j)  
   *      the access to these points is done by a LagrangeDOFHandler class.
   *    - a basis function phi_i of a neuman boundary point x_i 
   *      vanishes on the Robin-boundary and vice-versa. A basis function of an 
   *      interior point x_i vanishes on the boundary
   *    - The Dirichlet-Boundary is a closed set, in particular 
   *      point-evaluations indicate whether any of the adjacent edges or 
   *      faces are Dirichlet-boundaries
   *    - complete cell boundaries are of one type (except perhaps their 
   *      lower-codim boundaries): cog-evaluation indicates, whether 
   *      all Lagrange nodes on the intersection are Dirichlet-vertices. 
   *      And cog-evaluation indicates, 
   *      whether whole intersection is of one type, e.g. Neuman or Robin, for 
   *      integration over it. 
   *
   *  weak formulation of the above problem and restriction to the discrete 
   *  function with \f$ u_h := sum_j u_j phi_j \f$  leads to a system for the 
   *  DOF-vector u:
   *  
   *       \f$   M u = b \f$
   *
   *  with \f$
   *               
   *             /   kronecker(i,j)         if x_i is Dirichlet-LagrangePoint 
   *            /
   *    M_ij :=<       \int_\Omega     [a     grad(phi_j) ]^T  grad(phi_i) 
   *            \   -  \int_\Omega     [b     phi_j]^T         grad(phi_i)
   *             \  +  \int_\Omega     c          phi_i       phi_j
   *              \ +  \int_{\Gamma_R} alpha      phi_i       phi_j      otherwise
   *  
   *  and
   *
   *           /    g_D(x_i)               if x_i is Dirichlet-LagrangePoint
   *    b_i :=<   
   *           \      \int_\Omega     f   phi_i
   *            \   + \int_{\Gamma_N} g_N phi_i
   *             \  + \int_{\Gamma_R} g_R phi_i                        otherwise
   *
   *  \f$
   *  The right hand side is assumed to be assembled by another class, e.g.
   *  RhsAssembler, which is based on element-wise contributions
   *  by a RhsIntegrator class, etc.  
   *
   *  The matrix M has kronecker rows for all dirichlet points, but no 
   *  kronecker-columns. 
   *
   *  Optionally, the matrix and the right hand side can
   *  be processed to have also kronecker-columns, which is beneficial in case
   *  of a symmetric problem. This is done by calling the methods (the latter 
   *  possibly being repeated for changing or multiple rhsides.)
   *
   *       matrixKroneckerColumnsTreatment();
   *       rhsKroneckerColumnsTreatment(rhs);
   * 
   *  This results in
   *
   *     \f$           M_sym u = b_sym \f$
   *
   *  The new matrix has entries
   *
   *               /   kronecker(i,j)    if x_i OR x_j is Dirichlet-Lagr.Point 
   *   M_ij_sym :=<    M_ij              otherwise
   * 
   *  The new right hand side has entries:
   *
   *               /    b_i                  if x_i is Dirichlet-Lagr..Point
   *    b_i_sym :=<   
   *               \    b_i - sum_{x_j Dirichlet-Point} M_ij g_D(x_j)   otherwise
   *
   *  In this case, the deleted matrix entries are stored, as they are 
   *  required for every subsequent right hand side modification.
   *  This Kronecker-Column Treatment is performed by storing the following
   *  temporary objects:
   *
   *  \f$ d_dir := \f$ vector with 0 for non-Dirichlet DOFs, 1 for DirichletDOFs
   *  \f$ M_dir := \f$Null-matrix with all Dirichlet-Columns of M  - diag(d_dir)
   *             i.e. Dirichlet-Rows are completely zero
   *  \f$ g_D := \f$ vector with 0 for non-Dirichlet DOFs, 
               Dirichlet-Value for DirichletDOFs
   *
   *  Then, the symmetrization can compactly be written as
   *
   *   \f$ M_sym = M - M_dir \f$
   *   \f$ b_sym = b - M_dir * g_D \f$
   *
   *  The class depends on two template parameters, a SystemMatrixImp and an 
   *  ElementMatrixIntegratorImp class
   *
   *  The SystemMatrixImp class must provide a storage type for the global 
   *  operator matrix and some arithmetics and basic functionality. In 
   *  particular a print() method and an apply() method are assumed to exist.
   *  Currently the SparseRowMatrix satisfies the interface required for a
   *  SystemMatrixImp class.
   * 
   *  The ElementMatrixIntegratorImp class provides functionality for 
   *  computing a
   *  local element matrix on a given entity without dirichlet-treatment, i.e.
   *  by traversing the grid with the ElementMatrixImp class results in 
   *  the above matrix M_ij, with Dirichlet-values also integrated. After this
   *  Assembly, a Dirichlet-BndCorrection is performed.
   *
   *  Different operating modes are possible (currently however only 
   *  ASSEMBLED is implemented). In case of ASSEMBLED, the whole 
   *  global operator matrix is allocated and completely precomputed by 
   *  corresponding methods. In case of ON_THE_FLY, no complete matrix is
   *  allocated, but the matrix-vector multiplication is performed by 
   *  on-the-fly computation of the elementwise local matrices.
   */
  template< class SystemMatrixImp, class ElementMatrixIntegratorImp >
  class FEOp
  : public OEMSolver :: PreconditionInterface,
    public Operator< typename ElementMatrixIntegratorImp :: TraitsType :: DomainFieldType,
                     typename ElementMatrixIntegratorImp :: TraitsType :: RangeFieldType,
                     typename ElementMatrixIntegratorImp :: TraitsType :: DiscreteFunctionType,
                     typename ElementMatrixIntegratorImp :: TraitsType :: DiscreteFunctionType >
  {
  public:
    //! Type of matrix storage class used for global operator matrix
    typedef SystemMatrixImp SystemMatrixType;

    //! Type of ElementMatrixContributor
    typedef ElementMatrixIntegratorImp ElementMatrixIntegratorType;
  
    //! Operation mode: global allocation and multiplication or only 
    //! on-the-fly
    enum OpMode { ON_THE_FLY, ASSEMBLED };

    //! Mode of operation on the matrix/rhs concerning Dirichlet rows/cols
    // enum DirichletTreatmentMode { KRONECKER_ROWS, KRONECKER_ROWS_COLS };

    typedef typename ElementMatrixIntegratorType :: TraitsType TraitsType;
    typedef typename ElementMatrixIntegratorType :: ModelType ModelType;
    
    typedef typename TraitsType :: DiscreteFunctionType DiscreteFunctionType;
    typedef typename TraitsType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename TraitsType :: ElementMatrixType ElementMatrixType;
    typedef typename TraitsType :: IntersectionQuadratureType IntersectionQuadratureType;

  private: 
    //! reference to element matrix generator provided during initialization
    ElementMatrixIntegratorType &elementMatrixIntegrator_;

    //! the corresponding function_space 
    const DiscreteFunctionSpaceType &functionSpace_;
   
    //! pointer to the representing global matrix 
    mutable SystemMatrixType *matrix_; 

    //! flag indicating whether the global matrix is assembled 
    mutable bool matrix_assembled_;

    //! operator mode 
    const OpMode opMode_;

    //! maximal number of nonzeros per row in the global matrix. Used for 
    //! allocation. 
    const int maxNonZerosPerRow_;

    //! dirichlet-treatment mode 
    // DirichletTreatmentMode dirichletMode_;

    //! vector (matrix with single row), which indicates the Dirichlet-DOFs 
    mutable SparseRowMatrix< int > *isDirichletDOF_;

    //! flag indicating whether the DirichletDOF lookup table is assembled 
    mutable bool isDirichletDOF_assembled_;

    //! matrix for storage of the matrix columns, which are deleted 
    //! during symmetrization, but required for RHS modification 
    mutable SparseRowMatrix< double > *matrixDirichletColumns_;
   
    //! verbosity flag
    const int verbose_;

    const bool preconditionSSOR_;
 
    // //! pointers to storage of argument and destination, only required in 
    // LocalOperator Mode
    // const DiscreteFunctionType * arg_;
    // DiscreteFunctionType * dest_;

  public:
    /*!  constructor: Initialization of FEOp
     *
     *   Based on an existing instance of an elementmatrixintegrator 
     *   (which knows the model and the function space) the FEOp is initialized. 
     *   Operation mode must be selected as ASSEMBLED or ON_THE_FLY. The number of 
     *   nonzeros in the global matrix must be specified (only relevant in 
     *   ASSEMBLED-Mode). 
     * 
     *   \param elementMatrixIntegrator an instance of an element matrix integrator
     *
     *   \param opMode the operator mode 
     *
     *   \param maxNonZerosPerRow the maximum number of nonzeros per row in the 
     *          global matrix (default 50)
     *
     *   \param verbose set to zero
     *   \param preconditionSSOR parameter for precondition is set to false
     *   \return the initialized FEOp
     */
    FEOp( ElementMatrixIntegratorType &elementMatrixIntegrator,
          OpMode opMode = ASSEMBLED,
          // DirichletTreatmentMode dirichletMode = KRONECKER_ROWS,
          int maxNonZerosPerRow = 50,
          int verbose = 0,
          bool preconditionSSOR = false )
    : elementMatrixIntegrator_( elementMatrixIntegrator ),
      functionSpace_( elementMatrixIntegrator_.discreteFunctionSpace() ),
      matrix_( NULL ), 
      matrix_assembled_( false ),
      // arg_ ( NULL ), 
      // dest_ (NULL) , 
      opMode_( opMode ),
      // dirichletMode_( dirichletMode ),
      maxNonZerosPerRow_( maxNonZerosPerRow ),
      isDirichletDOF_( NULL ),
      isDirichletDOF_assembled_( false ),
      matrixDirichletColumns_( NULL ),
      verbose_( verbose ),
      preconditionSSOR_( preconditionSSOR )
    {
      if( verbose_ )
        std :: cout << "entered constructor of FEOp" << std :: endl;

      // class currently only implemented for non-symmetrized matrix/rhs!
      assert( opMode == ASSEMBLED );
    }
  
    /*! destructor: In case of allocation of global operator matrix
     *              it is deallocated, also the dirichletDOF lookup table.
     */
    ~FEOp () 
    {
      if( verbose_ )
          std :: cout << "entered destructor of FEOp" << std :: endl;

      if( matrix_ != NULL ) 
        delete matrix_;
      if( isDirichletDOF_ != NULL )
        delete isDirichletDOF_;
      if( matrixDirichletColumns_ != NULL ) 
        delete matrixDirichletColumns_;
    }

    /** \brief print matrix to output stream
     *
     *  \note This method only makes sense in ASSEMBLED mode.
     *
     *  \param[in]  out  stream to print matrix to (defaults to std::cout)
     */
    void print ( std :: ostream &out = std :: cout ) const 
    {
      if( verbose_ )
        std :: cout << "entered print() of FEOp" << std :: endl;
      
      assert( opMode_ == ASSEMBLED );
      systemMatrix().print( out );
    }

    /** \brief obtain assembled system matrix
     *
     *  \note This method only makes sense in ASSEMBLED mode.
     *
     *  \returns a reference to assembled system matrix
     */
    SystemMatrixType& systemMatrix () const
    {
      if( verbose_ )
        std :: cout << "entered systemMatrix() of FEOp" << std :: endl;
      
      assert( opMode_==ASSEMBLED );
      if( !matrix_assembled_ )
        assemble();
      return *matrix_;
    }
  
    /** \brief application operator
     *
     *   In case of ASSEMBLED-mode a global apply() on the matrix is called,
     *   i.e. a matrix-vector-multiplication with the vector arg is performed,
     *   the result stored in dest. This is an implicit requirement on the
     *   MatrixImp class!
     *   
     *   In case of ON_THE_FLY, the method is not yet implemented.
     *
     *   \param[in]   arg   reference on the argument for the matrix-vector
     *                      multiplication 
     *   \param[out]  dest  reference to storage for the destination vector
     */
    virtual void operator() ( const DiscreteFunctionType &arg,
                              DiscreteFunctionType &dest ) const
    {
      if( opMode_ == ASSEMBLED )
        systemMatrix().apply( arg, dest );
      else
        DUNE_THROW( NotImplemented,
                    "operator() in FEOP needs to be implemented for ON_THE_FLY!" );
    }

    /** \brief method required by oemsolvers */
    void finalize ()
    { 
    };

    /**  \brief perform grid-walkthrough and assemble global matrix
     * 
     *   If the matrix storage is not allocated, new storage is allocated by
     *   allocateSystemMatrix.
     *   The begin and end iterators are determined, 
     *   the assembling of the global matrix is initiated by call of 
     *   assembleOnGrid 
     *   and the Dirichlet-rows deleted by bndCorrectMatrixOnGrid. 
     *   The assemled flag is set. 
     *
     *   Method only makes sense in ASSEMBLED-mode
     */
    void assemble () const
    {
      if( verbose_ )
        std :: cout << "entered assemble() of FEOp" << std :: endl;

      assert( opMode_==ASSEMBLED );
        
      if( matrix_ == NULL )
        allocateSystemMatrix();

      matrix_->clear();
      // generate global matrix without dirichlet-treatment 
      assembleOnGrid();
      // generate kronecker-rows for dirichlet DOFs  
      bndCorrectMatrix();
    
      // in case of dirichletTreatment by Kronecker-kill: eliminate columns
      //if (dirichletMode_==KRONECKER_ROWS_COLS)
      //    matrixKroneckerColumnsTreatment();

      if( verbose_ )
        determineRealNonZeros();

      // in case of use of SSOR preconditioning, the matrix must not have
      // zero diagonal entries. Therefore zero-rows are transformed
      // to unit-rows and a global matrix diagonal check is performed.
      if( preconditionSSOR_ )
      {
        const int numRows = matrix_->rows();
        // set unit rows in empty rows and check diagonal for zeros
        for( int i = 0; i != numRows; ++i )
        {
          if( matrix_->numNonZeros( i ) == 0 )
            matrix_->set( i, i, 1.0 );
          assert( (*matrix_)( i, i ) != 0 );
        }
      }
       
      matrix_assembled_ = true;
    };

    /**  \brief mark the local variables to be no longer actual
     *
     *   method marks the global matrix and the Dirichlet-DOF-list for 
     *   reassembly. This must be called, if the underlying model or 
     *   element-matrix has changed after initialization of the FEOp. By this, 
     *   the next assemble() or apply() or operator() call of the FEOp will 
     *   invoke a new computation of these internal quantities
     */
    void markForReassembling ()
    {
      if( verbose_ )
        std :: cout << "entered markForReassembling() of FEOp" << std :: endl;
      isDirichletDOF_assembled_ = false;
      matrix_assembled_ = false;
    }
  
    /** \brief produce kronecker-columns in matrix
     *
     *   This method can be called after Matrix assembly. It changes the matrix
     *   M (given above), such that kronecker-columns are produced for
     *   Dirichlet DoFs. For symmetric problems this results in a symmetric
     *   matrix, which is beneficial for system-solvers. Note, that the
     *   right-hand side similarly has to be modified by
     *   rhsKroneckerColumnsTreatment in  order to give the identical result.
     *   The old matrix values are stored in a member variable in order to be
     *   able to process arbitrary many right hand side vectors afterwards.
     *
     *   The new matrix has entries
     *
     *               /   kronecker(i,j)     if x_i OR x_j is Dirichlet-Bnd.Point 
     *   M_ij_sym :=<    M_ij               otherwise
     */
    void matrixKroneckerColumnsTreatment() const
        {
          if (verbose_)
              std::cout << "entered matrixKroneckerColumnsTreatment of FEOp\n";
//          assert(isDirichletDOF_old_);
          assert(isDirichletDOF_);
          assert(matrix_);
          
          // delete storage if allocated
//          if (matrixDirichletColumns_old_)
//              delete(matrixDirichletColumns_old_);
          if (matrixDirichletColumns_)
              delete(matrixDirichletColumns_);

          if (verbose_>=2)
              std::cout << "deleted old storage for submatrix\n";
          
          // count new number of Dirichlet-DOFs
//          int nDirDOFs_old = 0;
//          for (int i=0; i!=matrix_->rows(); i++)
//              if ((*isDirichletDOF_old_)(0,i)) // so is Dirichlet-DOF
//                  nDirDOFs_old++;

          int nDirDOFs = 0;
          for (int i=0; i!=matrix_->rows(); i++)
              if ((*isDirichletDOF_)(0,i)) // so is Dirichlet-DOF
                  nDirDOFs++;
 
//          assert(nDirDOFs == nDirDOFs_old);
                   
          if (verbose_>=2)
              std::cout << "counted number of dirichletDOFs : " << 
                  nDirDOFs << " \n";

          // allocate new storage for to be deleted matrix entries
//          matrixDirichletColumns_old_ = new 
//                           OldSparseRowMatrix<double>
//              (matrix_->rows(), 
//               matrix_->cols(), 
//               maxNonZerosPerRow_);
// tooo much!!!:         nDirDOFs);

          matrixDirichletColumns_ = new 
                           SparseRowMatrix<double>
              (matrix_->rows(), 
               matrix_->cols(), 
               maxNonZerosPerRow_, 
               0.0);
          
          if (verbose_>=2)
              std::cout << "allocated new storage for submatrix\n";
          
//          assert(matrixDirichletColumns_old_);
          assert(matrixDirichletColumns_);
          
          // save matrix entries 
//          matrixDirichletColumns_old_->clear();
//          for (int i=0; i!=matrix_->rows(); i++)
//          {
//            if (i%100 == 0)
//                if (verbose_)
//                    std::cout << " processing row no " << i << " \n";
//            typename OldSparseRowMatrix<int>::ColumnIterator it = 
//                isDirichletDOF_old_->rbegin(0);
//            for (;it!=isDirichletDOF_old_->rend(0);++it)
//                matrixDirichletColumns_old_->set(i,it.col(),
//                                            (*matrix_)(i,it.col()));
//          }

          matrixDirichletColumns_->clear();
          
          for (int i=0; i!=matrix_->rows(); i++)
          {
            if (i%100 == 0)
                if (verbose_>=2)
                    std::cout << " processing row no " << i << " \n";
            // for each non-Dirichlet-row
            if (!(*isDirichletDOF_)(0,i)) 
            {
              // copy all dirichlet-DOF-entries from matrix              
              int numNonZeros = matrix_-> numNonZeros(i);
              for (int fakeCol =0; fakeCol!= numNonZeros; fakeCol++)
              {
                int realCol = matrix_->
                    realCol(i,fakeCol);
                if ((realCol != SparseRowMatrix<int>::defaultCol)
                    && ((*isDirichletDOF_)(0,realCol))) 
                    matrixDirichletColumns_->set(i, realCol,
                                                 (*matrix_)(i,realCol));
              }
            }
          }
           
          // delete matrix columns by subtracting matrices
                    
          int numNonZeros = isDirichletDOF_->numNonZeros(0);
          for (int fakeCol = 0; fakeCol!=numNonZeros; fakeCol++)
          {
            int realCol = isDirichletDOF_->realCol(0,fakeCol);
            if (realCol!=SparseRowMatrix<int>::defaultCol)
            { 
              if (verbose_>=2)
                  std::cout << " setting unitcolumn no " << 
                      realCol << " \n";
              matrix_->unitCol( realCol );    
            }
          }
          
          if (verbose_>=2)
              determineRealNonZeros();
        };
  
    /** \brief modify computed right-hand side for dirichlet nodes
     *
     *  This method is an optional postprocessing step. It changes the
     *  right-hand-side \f$b\f$ (assumed to be given as defined above), such
     *  that the symmetrized matrix \f$\tilde M\f$ (after application of
     *  matrixKroneckerColumnsTreatment) with the new rhs produces the same
     *  solution as the original matrix \f$M\f$ with the original rhs \f$b\f$.
     *  The new vector has the following entries:
     *  \f{displaymath}
     *  \tilde b_i = \left\lbrace\begin{array}{ll}
     *    b_i,
     *      & \enspace \mbox{if} \enspace x_i \in \Gamma_D,\\
     *    b_i - \sum_{x_j \in \Gamma_D} M_{ij} g_D( x_j )
     *      & \enspace \mbox{otherwise}.
     *  \end{array}\right.
     *  \f}
     *
     *  \note The values \f$g_D( x_j )\f$ are exactly assumed to be the values
     *  of the provided vector \f$b\f$. The boundary data is not reevaluated.
     *
     *  \param  rhs  a reference to the assembled right-hand side \f$b\f$
     */
    void rhsKroneckerColumnsTreatment(DiscreteFunctionType& rhs) const
        {
          assert(isDirichletDOF_);
          assert(matrixDirichletColumns_);
          
          assert ( rhs.size() == isDirichletDOF_->cols() );
          
          if (verbose_)
              std::cout << "entered rhsKroneckerColumnsTreatment\n";

//          saveSparseMatrixBinary("dircol_matrix.bin",*matrixDirichletColumns_);
//          saveDofVectorBinary("rhs_before_symm.bin",rhs);

          // allocate storage for update vector
          
          const typename DiscreteFunctionType::DiscreteFunctionSpaceType & 
              discFuncSpace = rhs.space();
          
          DiscreteFunctionType update = 
              DiscreteFunctionType( "update", discFuncSpace);

//          saveDofVectorBinary("update_before_mult.bin",rhs);
          
          // determine update = M_dir * b = M_dir * 
          matrixDirichletColumns_->apply(rhs, update );

//          saveDofVectorBinary("update_after_mult.bin",rhs);
          
          // compute rhs = b_sym = b - M_dir * b = rhs - update
//          saveDofVectorBinary("rhs_before_subtracting.bin",rhs);
          rhs -=update;
//          saveDofVectorBinary("rhs_after_subtracting.bin",rhs);
          
          if (verbose_>=2)
              std::cout << " finished\n";
          
//          std::cout << "finished Rhs-KroneckerColumnTreatment\n";
        };



    //! method required for preconditioning in bicgstab, to be provided for 
    //! PreconditioningInterface
    bool hasPreconditionMatrix() const 
    { 
      return preconditionSSOR_; 
    }

    //! method required for preconditioning in bicgstab, to be provided for
    //! PreconditioningInterface
    const SystemMatrixType& preconditionMatrix () const
    {
      return systemMatrix();
    }
  
  private:
    /*!  allocateSystemMatrix: allocation of a new system matrix
     *
     *   \note SystemMatrixType requires a constructor with the following syntax:
     *   SystemMatrixType( nrows, ncols, nonzeros_per_row ).
     *   
     *   \note Currently the nonzeros per row must be specified in the
     *   constructor of FEOp. This could be improved, e.g. by a Traitsclass as
     *   template-parameter, which  contains the SystemMatrixType and maxnumbernonzeros
     */
    void allocateSystemMatrix () const  
    { 
      if( verbose_ )
        std :: cout << "entered allocateSystemMatrix() of FEOp" << std :: endl;

      const size_t size = functionSpace_.size();
      if( verbose_ )
        std :: cout << "allocating matrix of size " << size
                    << " times " << maxNonZerosPerRow_
                    << " nonzeros" << std :: endl;
          
      matrix_ = new SystemMatrixType( size, size, maxNonZerosPerRow_ );
      assert( matrix_ != NULL );
    }
    
    /*! 
     *   assembleOnGrid: perform grid walkthrough and assemble matrix
     *
     *   For each element, the local element matrix is determined into the 
     *   given local matrix storage and distributed into the global matrix.
     *   Distribution is performed by an add(row,col,val) method on the 
     *   global matrix class.
     */
    void assembleOnGrid () const
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
        BaseFunctionSetType;

      if( verbose_ )
        std::cout << "entered assembleOnGrid() of FEOp" << std :: endl;
      
      ElementMatrixType elementMatrix;
      // run through grid and add up local contributions
      const IteratorType endit = functionSpace_.end();
      for( IteratorType it = functionSpace_.begin(); it != endit; ++it )
      {
        const BaseFunctionSetType baseFunctionSet
          = functionSpace_.baseFunctionSet( *it );
        const int numBaseFunctions = baseFunctionSet.numBaseFunctions();
            
        // setup local element matrix 
        // (size check is performed in elementmatrix construction)
        elementMatrix.clear();
        elementMatrixIntegrator_.addElementMatrix( *it, elementMatrix, 1 );
            
        for( int i = 0; i < numBaseFunctions; ++i ) 
        { 
          const int row = functionSpace_.mapToGlobal( *it , i );
          for( int j = 0; j < numBaseFunctions; ++j ) 
          {
            const int col = functionSpace_.mapToGlobal( *it , j );    
            matrix_->add( row, col, elementMatrix( i, j ) );
          }
        }
      }
    }

    /*! searchDirichletDOFs: determine lookup table for dirichlet-values
     *
     *   method fills the local vector isDirichletDOF_old_ as multiple operations with
     *   Dirichlet-boundaries are necessary, e.g. matrix-boundary treatment, 
     *   NOT rhs-assembly but later symmetrization of the system. By this multiple 
     *   global grid walkthroughs can be prevented and replaced by single run 
     *   over the lookup-table.
     */
    void searchDirichletDOFs () const
    {
      if( verbose_ >= 1 )
        std :: cout << "entered searchDirichletDOFs() of FEOp" << std :: endl;
      
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType :: GridType GridType; 
      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename GridPartType :: IntersectionIteratorType 
        IntersectionIteratorType;
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      // codimension of faces
      enum { faceCodim = 1 };
      
      // type of Lagrange point set
      // note: DiscreteFunctionSpaceType must be a Lagrange discrete function
      //       space!
      typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
        LagrangePointSetType;
      // type of iterator over all DoFs of a face
      typedef typename LagrangePointSetType :: template Codim< faceCodim >
                                            :: SubEntityIteratorType
        FaceDofIteratorType;

      // allocate isDirichlet-Vector if not already done
      if( this->isDirichletDOF_ == NULL ) {
        const unsigned int size = functionSpace_.size();
        isDirichletDOF_  = new SparseRowMatrix< int >( 1, size, size, 0 );
      }
      assert( this->isDirichletDOF_ != NULL );
      if( verbose_ >= 1 )
        std :: cout << "allocated isDirichletDOF." << std :: endl;
          
      const GridPartType &gridPart = functionSpace_.gridPart();

      if( verbose_ >= 1 )
        std :: cout << "starting loop over grid." << std :: endl;
 
      IteratorType it = functionSpace_.begin();
      const IteratorType endit = functionSpace_.end(); 
      for( int numit = 0; it != endit; ++it, ++numit ) {
        if( (verbose_ >= 2) && (numit % 100 == 0) )
          std :: cout << "processing entity no " << numit << std :: endl;

        const EntityType &entity = *it;

        IntersectionIteratorType nit = gridPart.ibegin( entity );
        const IntersectionIteratorType endnit = gridPart.iend( entity );
        for( ; nit != endnit; ++nit ) {
          // make sure we are on the boundary
          if( !nit.boundary() )
            continue;
              
          if( elementMatrixIntegrator_.model().boundaryType( nit ) != ModelType :: Dirichlet )
            continue;

          const int faceNumber = nit.numberInSelf();
          
          const LagrangePointSetType &lagrangePointSet
            = functionSpace_.lagrangePointSet( entity );
          
          FaceDofIteratorType faceIt
            = lagrangePointSet.template beginSubEntity< faceCodim >( faceNumber );
          const FaceDofIteratorType faceEndIt
            = lagrangePointSet.template endSubEntity< faceCodim >( faceNumber );
          for( ; faceIt != faceEndIt; ++faceIt ) {
            const int row = functionSpace_.mapToGlobal( entity, *faceIt );
            isDirichletDOF_->set( 0, row, 1 );
          }
        }
      }
      
      isDirichletDOF_assembled_ = true;
        
      if( verbose_>=2 )
        std :: cout << "finished searchDirichletDOFs" << std :: endl;
    }
  
    /*! 
     *   bndCorrectMatrix: treatment of Dirichlet-DOFS
     *
     *   delete rows for dirichlet DOFS, setting diagonal 
     *   element to 1. This is reasonable, as Lagrange Basis is implicitly 
     *   assumed, 
     *   so the RHS being the exact dirichlet-values, therefore, the matrix 
     *   row must be a unit-row.
     */
    void bndCorrectMatrix () const
    {
      if( verbose_ )
        std :: cout << "entered bndCorrectMatrix() of FEOp" << std :: endl;

      if( !isDirichletDOF_assembled_ )
        searchDirichletDOFs();
          
      // eliminate the Dirichlet rows by converting to unit-rows      
      const int numNonZeros = isDirichletDOF_->numNonZeros( 0 );
      for( int fakeCol = 0; fakeCol < numNonZeros; ++fakeCol )
      {
        const int realCol = isDirichletDOF_->realCol( 0, fakeCol );
        if( realCol != SparseRowMatrix< int > :: defaultCol )
        {
          if( (*isDirichletDOF_)( 0, realCol ) )
          {
            if( verbose_ >= 2 )
              std :: cout << " setting unit row " << realCol << std :: endl;
            matrix_->unitRow( realCol );
          }
        } 
      }
          
      if( verbose_ >= 2 )
        std :: cout << " end of bndCorrectMatrix" << std :: endl;
    }

  /** \brief obtain the maximum number of nonzero entries in the system matrix
   *
   *  \returns the maximum real number of nonzero entries in the system matrix
   */
  int determineRealNonZeros() const
  {
    // search number of nonzeros
    if (verbose_)
        std::cout << "entered determineRealNonZeros() of FEOp \n";  

    int maxnonzeros = -1;
    
    for (int i=0; i!=matrix_->rows(); i++)
    {
      if ((verbose_>=2) && (i%100==0)) 
          std::cout << " counting nonzeros in row " << i <<" \n";  
      
      int nonzeros = 0;
      int numNonZeros = matrix_->numNonZeros(i);
      
      for (int fakeCol=0; fakeCol!=numNonZeros ; fakeCol++ )
      {
        int realCol = matrix_->realCol(i,fakeCol);
        if (realCol != SystemMatrixType::defaultCol)
            if ((*matrix_)(i,realCol)!=0.0)
                nonzeros++;
      }
      
//      typename SystemMatrixType::ColumnIterator 
//          it = matrix_->rbegin(i); 
//      typename SystemMatrixType::ColumnIterator 
//          endit = matrix_->rend(i);
      
//      for (;it!=endit;++it)
//          if (*it!=0.0) nonzeros++;
      
//      for (int j=0; j!=matrix_->cols(); j++)
//          if ((*matrix_)(i,j)!=0.0) nonzeros++;
//      for (int j=0; j!=matrix_->cols(); j++)
//          if ((*matrix_)(i,j)!=0.0) nonzeros++;
      
      if (nonzeros > maxnonzeros)
          maxnonzeros = nonzeros;
    }
    std::cout << " current real nonzeros per row:" << maxnonzeros << " \n";
    
    return maxnonzeros;
  }
}; // end class FEOp

} // end namespace


#endif
