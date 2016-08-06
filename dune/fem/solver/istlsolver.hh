#ifndef DUNE_FEM_ISTLSOLVERS_HH
#define DUNE_FEM_ISTLSOLVERS_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>


#if HAVE_DUNE_ISTL
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/fem/operator/linear/istladapter.hh>
#include <dune/fem/operator/linear/istloperator.hh>

namespace Dune
{

  namespace Fem
  {

    //=====================================================================
    // Implementation for ISTL-matrix based operator
    //=====================================================================

    /** @ingroup OEMSolver
        @{
    **/

    // DefaultSolverCaller
    // -------------------

    template< class OperatorImp, class DiscreteFunction, class SolverCaller >
    struct DefaultSolverCaller
    {
      template <class MatrixAdapter>
      static std::pair< int, double >
      call ( const OperatorImp &op,
             MatrixAdapter& matrix,
             const DiscreteFunction &arg, DiscreteFunction &dest,
             double reduction, double absLimit, int maxIter, bool verbose,
             const ParameterReader &parameter )
      {
        typedef typename DiscreteFunction :: DofStorageType BlockVectorType;

        // verbose only in verbose mode and for rank 0
        int verb = (verbose && (dest.space().gridPart().comm().rank() == 0)) ? 2 : 0;

        if( absLimit < std::numeric_limits< double >::max() )
        {
          const double residuum = matrix.residuum( arg.blockVector(), dest.blockVector() );
          reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

          if( verbose && (dest.space().gridPart().comm().rank() == 0) )
            std::cout << SolverCaller::name() <<": reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit << std::endl;
        }

        typename SolverCaller :: SolverType solver( matrix, dest.scalarProduct(),
                                                    matrix.preconditionAdapter(), reduction, maxIter, verb );

        // copy right hand side since ISTL is overwriting it
        BlockVectorType rhs( arg.blockVector() );

        InverseOperatorResult returnInfo;

        // call solver
        solver.apply( dest.blockVector(), rhs, returnInfo );

        return std::pair< int, double > ( returnInfo.iterations, matrix.averageCommTime() );
      }

    };


    // ISTLInverseOp
    // Baseclass for each ISTL solver
    // -------------

    template< class DF, class Op, class SolverCaller >
    struct ISTLInverseOp
    : public Operator< DF, DF >
    {
    public:
      typedef DF DiscreteFunctionType;
      typedef DiscreteFunctionType  DestinationType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef Op OperatorType;
      typedef SolverCaller SolverCallerType;

      typedef Dune::Fem::ISTLMatrixFreeOperatorAdapter< OperatorType >                      ISTLMatrixFreeAdapterType;
      typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType >   AssembledOperatorType;

      /** \brief constructor
       *
       *  \param[in] op Mapping describing operator to invert
       *  \param[in] reduction reduction epsilon
       *  \param[in] absLimit absolut limit of residual (not used here)
       *  \param[in] maxIter maximal iteration steps
       *  \param[in] verbose verbosity
       *
       *  \note ISTL BiCG-stab only uses the relative reduction.
       */
      ISTLInverseOp ( const OperatorType &op,
                      double reduction, double absLimit, int maxIter, bool verbose,
                      const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        matrixOp_( dynamic_cast<const AssembledOperatorType*> (&op_) ),
        reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 ),
        averageCommTime_( 0.0 ),
        parameter_( parameter )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      ISTLInverseOp ( const OperatorType &op,
                      double reduction, double absLimit, int maxIter,
                      const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        matrixOp_( dynamic_cast<const AssembledOperatorType*> (&op_) ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        averageCommTime_( 0.0 ),
        parameter_( parameter )
      {}

      ISTLInverseOp ( const OperatorType &op,
                      double reduction, double absLimit,
                      const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        matrixOp_( dynamic_cast<const AssembledOperatorType*> (&op_) ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( std::numeric_limits< int >::max() ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        averageCommTime_( 0.0 ),
        parameter_( parameter )
      {}

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {}

      void finalize () const
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver:"<< SolverCallerType::name() << ",  eps = " << reduction_ ;
        out  << "\\\\ \n";
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        // if op_ is an instance of ISTLLinearOperator (i.e. assembled) we can
        // use the corresponding matrix adapter, otherwise we assume it's a
        // matrix free implementation

        std::pair< int, double > info;
        if( matrixOp_ )
        {
          typedef typename AssembledOperatorType :: BaseType  MatrixObjectType;
          const MatrixObjectType& matrixObj = matrixOp_->systemMatrix() ;

          typedef ISTLMatrixAdapterFactory< MatrixObjectType > ISTLMatrixAdapterFactoryType;
          auto matrixAdapterPtr = ISTLMatrixAdapterFactoryType :: matrixAdapter( matrixObj );
          info = SolverCallerType::call( op_, *matrixAdapterPtr,
                                         arg, dest, reduction_, absLimit_, maxIter_, verbose_, parameter_ );
        }
        else
        {
          ISTLMatrixFreeAdapterType matrixAdapter( op_, arg.space(), dest.space() );
          info = SolverCallerType::call( op_, matrixAdapter,
                                         arg, dest, reduction_, absLimit_, maxIter_, verbose_, parameter_ );
        }

        iterations_ = info.first;
        averageCommTime_ = info.second;
      }

      // return number of iterations
      int iterations() const
      {
        return iterations_;
      }

      //! return accumulated communication time
      double averageCommTime() const
      {
        return averageCommTime_;
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void operator() ( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply( arg,dest );
      }

    private:
      const OperatorType &op_;
      const AssembledOperatorType* matrixOp_;
      double reduction_;
      double absLimit_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;
      mutable double averageCommTime_;
      ParameterReader parameter_;
    };


    // LoopSolverCaller
    // --------------

    template< class OperatorImp, class DiscreteFunction >
    struct LoopSolverCaller
    : public DefaultSolverCaller< OperatorImp, DiscreteFunction, LoopSolverCaller< OperatorImp, DiscreteFunction > >
    {
      typedef LoopSolver< typename DiscreteFunction :: DofStorageType > SolverType;

      static inline std::string name ()
      {
        std::string name = "ISTL LoopSolver";
        return name;
      }
    };


    // ISTLLoopOp
    // --------------

    /** \brief LoopSolver scheme for block matrices (BCRSMatrix)
               and block vectors (BVector) from DUNE-ISTL.
        \note Op defaults to Dune::Fem::Operator.
     */
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    struct ISTLLoopOp
    : public ISTLInverseOp< DF, Op, LoopSolverCaller< Op, DF > >
    {
      typedef ISTLInverseOp< DF, Op, LoopSolverCaller< Op, DF > > BaseType;
    public:

      template< class ... Args >
      ISTLLoopOp ( Args && ... args ) : BaseType( std::forward< Args >( args ) ... )  {}
    };


    // MINResSolverCaller
    // ------------------

    template< class OperatorImp, class DiscreteFunction >
    struct MINResSolverCaller
    : public DefaultSolverCaller< OperatorImp, DiscreteFunction, MINResSolverCaller< OperatorImp, DiscreteFunction > >
    {
      typedef MINRESSolver< typename DiscreteFunction :: DofStorageType > SolverType;

      static inline std::string name ()
      {
        std::string name = "ISTL MINRESSolver";
        return name;
      }
    };


    // ISTLMINResOp
    // --------------

    /** \brief MINRes scheme for block matrices (BCRSMatrix)
        and block vectors (BVector) from DUNE-ISTL.
        \note Op defaults to Dune::Fem::Operator.
     */
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    struct ISTLMINResOp
    : public ISTLInverseOp< DF, Op, MINResSolverCaller< Op, DF > >
    {
      typedef ISTLInverseOp< DF, Op, MINResSolverCaller< Op, DF > > BaseType;
    public:
      template< class ... Args >
      ISTLMINResOp ( Args && ... args ) : BaseType( std::forward< Args >( args ) ... )  {}
    };


    // BiCGSTABSolverCaller
    // --------------------

    template< class OperatorImp, class DiscreteFunction >
    struct BiCGSTABSolverCaller
    : public DefaultSolverCaller< OperatorImp, DiscreteFunction, BiCGSTABSolverCaller< OperatorImp, DiscreteFunction > >
    {
      typedef BiCGSTABSolver< typename DiscreteFunction :: DofStorageType > SolverType;

      static inline std::string name ()
      {
        std::string name = "ISTL BiCGSTABSolver";
        return name;
      }
    };


    // ISTLBICGSTABOp
    // --------------

    /** \brief BICG-stab scheme for block matrices (BCRSMatrix)
        and block vectors (BVector) from DUNE-ISTL.
        \note Op defaults to Dune::Fem::Operator.
     */
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    struct ISTLBICGSTABOp
    : public ISTLInverseOp< DF, Op, BiCGSTABSolverCaller< Op, DF > >
    {
      typedef ISTLInverseOp< DF, Op, BiCGSTABSolverCaller< Op, DF > > BaseType;
    public:
      template< class ... Args >
      ISTLBICGSTABOp ( Args && ... args ) : BaseType( std::forward< Args >( args ) ... )  {}
    };


    // GMResSolverCaller
    // -----------------

    template< class OperatorImp, class DiscreteFunction >
    struct GMResSolverCaller
    {
      template <class MatrixAdapter>
      static std::pair< int, double >
      call ( const OperatorImp &op,
             MatrixAdapter matrix,
             const DiscreteFunction &arg, DiscreteFunction &dest,
             double reduction, double absLimit, int maxIter, bool verbose,
             const ParameterReader &parameter )
      {
        int restart = parameter.getValue< int >( "istl.gmres.restart", 5 );

        typedef typename DiscreteFunction :: DofStorageType BlockVectorType;

        // verbose only in verbose mode and for rank 0
        int verb = (verbose && (dest.space().gridPart().comm().rank() == 0)) ? 2 : 0;

        if( absLimit < std::numeric_limits< double >::max() )
        {
          const double residuum = matrix.residuum( arg.blockVector(), dest.blockVector() );
          reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

          if( verbose && (dest.space().gridPart().comm().rank() == 0) )
            std::cout << "ISTL GMRes-Solver: reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit << std::endl;
        }

        RestartedGMResSolver< BlockVectorType >
          solver( matrix, matrix.scp(), matrix.preconditionAdapter(), reduction, restart, maxIter, verb );

        // copy right hand side since ISTL is overwriting it
        BlockVectorType rhs( arg.blockVector() );

        InverseOperatorResult returnInfo;

        // call solver
        solver.apply( dest.blockVector(), rhs, returnInfo );

        return std::pair< int, double > ( returnInfo.iterations, matrix.averageCommTime() );
      }

      static std::string name ()
      {
        std::string name = "ISTL GMRes-Solver";
        return name;
      }
    };


    // ISTLGMResOp
    // -----------

    /** \brief GMRes scheme for block matrices (BCRSMatrix)
     *         and block vectors (BVector) from dune-istl
     *  \note Op defaults to Dune::Fem::Operator.
     */
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    struct ISTLGMResOp
    : public ISTLInverseOp< DF, Op, GMResSolverCaller< Op, DF > >
    {
      typedef ISTLInverseOp< DF, Op, GMResSolverCaller< Op, DF > > BaseType;
    public:
      template< class ... Args >
      ISTLGMResOp ( Args && ... args ) : BaseType( std::forward< Args >( args ) ... )  {}
    };


    // CGSolverCaller
    // --------------

    template< class OperatorImp, class DiscreteFunction >
    struct CGSolverCaller
    : public DefaultSolverCaller< OperatorImp, DiscreteFunction, CGSolverCaller< OperatorImp, DiscreteFunction > >
    {
      typedef CGSolver< typename DiscreteFunction :: DofStorageType > SolverType;

      static inline std::string name ()
      {
        std::string name = "ISTL CGSolver";
        return name;
      }
    };


    // ISTLCGOp
    // --------

    /** \brief BICG-stab scheme for block matrices (BCRSMatrix)
               and block vectors (BVector) from dune-istl.
        \note Op defaults to Dune::Fem::Operator.
     */
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    struct ISTLCGOp
    : public ISTLInverseOp< DF, Op, CGSolverCaller< Op, DF > >
    {
      typedef ISTLInverseOp< DF, Op, CGSolverCaller< Op, DF > > BaseType;
    public:
      template< class ... Args >
      ISTLCGOp ( Args && ... args ) : BaseType( std::forward< Args >( args ) ... )  {}
    };


    // SuperLUSolverCaller
    // -------------------

    template< class OperatorImp, class DiscreteFunction >
    struct SuperLUSolverCaller
    {
      template <class MatrixAdapterDummy>
      static std::pair< int, double >
      call ( const OperatorImp &op,
             MatrixAdapterDummy dummy,
             const DiscreteFunction &arg, DiscreteFunction &dest,
             double reduction, double absLimit, int maxIter, bool verbose,
             const ParameterReader &parameter )
      {
        typedef typename OperatorImp :: MatrixAdapterType MatrixAdapterType;
        MatrixAdapterType matrix  = op.systemMatrix().matrixAdapter();

        static_assert( std::is_same< MatrixAdapterDummy, MatrixAdapterType > :: value,
                       "SuperLU only works with assembled operators" );

        InverseOperatorResult returnInfo;
#if HAVE_SUPERLU
        // create solver
        typedef typename MatrixAdapterType :: MatrixType MatrixType;
        typedef typename MatrixType :: BaseType ISTLMatrixType;
        SuperLU< ISTLMatrixType > solver( matrix.getmat(), verbose );

        typedef typename DiscreteFunction :: DofStorageType BlockVectorType;
        BlockVectorType rhs( arg.blockVector() );
        // solve the system
        solver.apply( dest.blockVector(), rhs, returnInfo );
#else
        DUNE_THROW(NotImplemented,"SuperLU solver not found in configure, please re-configure");
#endif

        // get information
        std::pair< int, double > p( returnInfo.iterations, matrix.averageCommTime() );
        return p;
      }

      static std::string name ()
      {
        std::string name = "ISTL SuperLU";
        return name;
      }
    };


    // ISTLSuperLUOp
    // --------------

    /** \brief SuperLU solver for block matrices (BCRSMatrix)
        and block vectors (BVector) from DUNE-ISTL.

        The solver is based in the
        well known <a href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU
        package</a>.
    */
    template< class DF, class Op >
    struct ISTLSuperLU
    : public ISTLInverseOp< DF, Op, SuperLUSolverCaller< Op, DF > >
    {
      typedef ISTLInverseOp< DF, Op, SuperLUSolverCaller< Op, DF > > BaseType;
    public:
      template< class ... Args >
      ISTLSuperLU ( Args && ... args ) : BaseType( std::forward< Args >( args ) ... )  {}
    };


  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLSOLVERS_HH
