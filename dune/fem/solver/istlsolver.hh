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

    struct ISTLSolverCall
    {
      template <class SolverType, class BlockVectorType>
      static int  
      solve( SolverType& solver, 
             const BlockVectorType& arg, 
             BlockVectorType& dest ) 
      {
        // copy right hand side since ISTL is overwriting it 
        BlockVectorType rhs( arg );

        InverseOperatorResult returnInfo;

        // call solver 
        solver.apply( dest, rhs, returnInfo );

        // get information 
        return returnInfo.iterations; 
      }
    };


    // ISTLMINResOp
    // --------------

    /** \brief MINRES scheme for block matrices (BCRSMatrix) 
        and block vectors (BVector) from DUNE-ISTL. */
    template< class DF, class Op >
    struct ISTLMINResOp
    : public Operator< DF, DF >
    {
    public:  
      typedef DF DiscreteFunctionType;
      typedef DiscreteFunctionType  DestinationType;
      typedef Op OperatorType;

    private:
      template <class OperatorImp, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static std::pair< int, double >
        call ( const OperatorImp &op,
               const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
               double reduction, double absLimit, int maxIter, bool verbose )
        {
          return solve( op.systemMatrix(), arg, dest, reduction, absLimit, maxIter, verbose );
        }

        template< class MatrixObjType, class DiscreteFunctionImp >
        static std::pair< int, double >
        solve ( const MatrixObjType &mObj,
                const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
                double reduction, double absLimit, int maxIter, bool verbose )
        {
          typedef typename MatrixObjType::MatrixAdapterType MatrixAdapterType;
          MatrixAdapterType matrix = mObj.matrixAdapter();
          
          typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;

          // verbose only in verbose mode and for rank 0 
          int verb = (verbose && (dest.space().gridPart().comm().rank() == 0)) ? 2 : 0;
           
          if( 0 < absLimit && absLimit < std::numeric_limits< double >::max() )
          {
            const double residuum = matrix.residuum( arg.blockVector(), dest.blockVector() );
            reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

            if( verbose ) 
              std::cout << "ISTL MINRes-Solver: reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit << std::endl;
          }

          MINRESSolver< BlockVectorType >
            solver( matrix, matrix.scp(), matrix.preconditionAdapter(), reduction, maxIter, verb );

          // call solver and return info
          int iter = ISTLSolverCall::solve( solver, arg.blockVector(), dest.blockVector() );
          return std::pair< int, double > ( iter, matrix.averageCommTime() );
        }
      };

    public:
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
      ISTLMINResOp ( const OperatorType &op,
                     double  reduction,
                     double absLimit,
                     int maxIter,
                     bool verbose )
      : op_( op ),
        reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      ISTLMINResOp ( const OperatorType &op,
                     double reduction,
                     double absLimit,
                     int maxIter = std::numeric_limits< int >::max() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {}

      void finalize () const
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: ISTL MINRes,  eps = " << reduction_ ;
        out  << "\\\\ \n";
      }

      /** \brief solve the system 
          \param[in] arg right hand side 
          \param[out] dest solution 
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        typedef SolverCaller< OperatorType, true > Caller;

        std::pair< int, double > info
          = Caller::call( op_, arg, dest, reduction_, absLimit_, maxIter_, verbose_ );

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
        apply(arg,dest);
      }

    private:
      const OperatorType &op_;
      double reduction_;
      double absLimit_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;
      mutable double averageCommTime_;
    }; 


    // ISTLBICGSTABOp
    // --------------

    /** \brief BICG-stab scheme for block matrices (BCRSMatrix) 
        and block vectors (BVector) from DUNE-ISTL. */
    template< class DF, class Op >
    struct ISTLBICGSTABOp
    : public Operator< DF, DF >
    {
    public:  
      typedef DF DiscreteFunctionType;
      typedef DiscreteFunctionType  DestinationType;
      typedef Op OperatorType;

    private:
      template <class OperatorImp, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static std::pair< int, double >
        call ( const OperatorImp &op,
               const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
               double reduction, double absLimit, int maxIter, bool verbose )
        {
          return solve( op.systemMatrix(), arg, dest, reduction, absLimit, maxIter, verbose );
        }

        template< class MatrixObjType, class DiscreteFunctionImp >
        static std::pair< int, double >
        solve ( const MatrixObjType &mObj,
                const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
                double reduction, double absLimit, int maxIter, bool verbose )
        {
          typedef typename MatrixObjType::MatrixAdapterType MatrixAdapterType;
          MatrixAdapterType matrix = mObj.matrixAdapter();
          
          typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;

          // verbose only in verbose mode and for rank 0 
          int verb = (verbose && (dest.space().gridPart().comm().rank() == 0)) ? 2 : 0;
           
          if( 0 < absLimit && absLimit < std::numeric_limits< double >::max() )
          {
            const double residuum = matrix.residuum( arg.blockVector(), dest.blockVector() );
            reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

            if( verbose ) 
              std::cout << "ISTL BiCGSTAB-Solver: reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit << std::endl;
          }

          BiCGSTABSolver< BlockVectorType >
            solver( matrix, matrix.scp(), matrix.preconditionAdapter(), reduction, maxIter, verb );

          // call solver and return info
          int iter = ISTLSolverCall::solve( solver, arg.blockVector(), dest.blockVector() );
          return std::pair< int, double > ( iter, matrix.averageCommTime() );
        }
      };

    public:
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
      ISTLBICGSTABOp ( const OperatorType &op,
                       double  reduction,
                       double absLimit,
                       int maxIter,
                       bool verbose )
      : op_( op ),
        reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      ISTLBICGSTABOp ( const OperatorType &op,
                       double reduction,
                       double absLimit,
                       int maxIter = std::numeric_limits< int >::max() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {}

      void finalize () const
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: ISTL BiCG-STAB,  eps = " << reduction_ ;
        out  << "\\\\ \n";
      }

      /** \brief solve the system 
          \param[in] arg right hand side 
          \param[out] dest solution 
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        typedef SolverCaller< OperatorType, true > Caller;

        std::pair< int, double > info
          = Caller::call( op_, arg, dest, reduction_, absLimit_, maxIter_, verbose_ );

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
        apply(arg,dest);
      }

    private:
      const OperatorType &op_;
      double reduction_;
      double absLimit_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;
      mutable double averageCommTime_;
    }; 



    // ISTLGMResOp
    // -----------

    /** \brief GMRes scheme for block matrices (BCRSMatrix)
     *         and block vectors (BVector) from dune-istl
     */
    template< class DF, class Op >
    struct ISTLGMResOp
    : public Operator< DF, DF >
    {
    public:  
      typedef DF DiscreteFunctionType;
      typedef Op OperatorType;
      typedef DiscreteFunctionType  DestinationType;

    private:
      template <class OperatorImp, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static std::pair< int, double >
        call ( const OperatorImp &op,
               const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
               double reduction, double absLimit, int restart, int maxIter, bool verbose )
        {
          return solve( op.systemMatrix(), arg, dest, reduction, absLimit, restart, maxIter, verbose );
        }

        template< class MatrixObjType, class DiscreteFunctionImp >
        static std::pair< int, double >
        solve ( const MatrixObjType &mObj,
                const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
                double reduction, double absLimit, int restart, int maxIter, bool verbose )
        {
          typedef typename MatrixObjType::MatrixAdapterType MatrixAdapterType;
          MatrixAdapterType matrix = mObj.matrixAdapter();
          
          typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;

          // verbose only in verbose mode and for rank 0 
          int verb = (verbose && (dest.space().gridPart().comm().rank() == 0)) ? 2 : 0;
           
          if( absLimit < std::numeric_limits< double >::max() )
          {
            const double residuum = matrix.residuum( arg.blockVector(), dest.blockVector() );
            reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

            if( verbose ) 
              std::cout << "ISTL GMRes-Solver: reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit << std::endl;
          }

          RestartedGMResSolver< BlockVectorType >
            solver( matrix, matrix.scp(), matrix.preconditionAdapter(), reduction, restart, maxIter, verb );

          // call solver and return info
          int iter = ISTLSolverCall::solve( solver, arg.blockVector(), dest.blockVector() );
          return std::pair< int, double > ( iter, matrix.averageCommTime() );
        }
      };

    public:
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
      ISTLGMResOp ( const OperatorType &op,
                    double  reduction,
                    double absLimit,
                    int maxIter,
                    bool verbose )
      : op_( op ),
        reduction_( reduction ),
        absLimit_( absLimit ),
        restart_( Parameter::getValue< int >( "istl.gmres.restart", 5 ) ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      ISTLGMResOp ( const OperatorType &op,
                    double reduction,
                    double absLimit,
                    int maxIter = std::numeric_limits< int >::max() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        restart_( Parameter::getValue< int >( "istl.gmres.restart", 5 ) ),
        maxIter_( maxIter ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {}

      void finalize () const
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: ISTL GMRes,  eps = " << reduction_ ;
        out  << "\\\\ \n";
      }

      /** \brief solve the system 
          \param[in] arg right hand side 
          \param[out] dest solution 
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        typedef SolverCaller< OperatorType, true > Caller;

        std::pair< int, double > info
          = Caller::call( op_, arg, dest, reduction_, absLimit_, restart_, maxIter_, verbose_ );

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
        apply(arg,dest);
      }

    private:
      const OperatorType &op_;
      double reduction_;
      double absLimit_;
      int restart_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;
      mutable double averageCommTime_;
    }; 



    // ISTLCGOp
    // --------

    /** \brief BICG-stab scheme for block matrices (BCRSMatrix) 
    and block vectors (BVector) from dune-istl. */
    template< class DF, class Op >
    struct ISTLCGOp
    : public Operator< DF, DF >
    {
    public:  
      typedef DF DiscreteFunctionType;
      typedef Op OperatorType;
      typedef DiscreteFunctionType  DestinationType;

    private:
      template< class OperatorImp, bool hasPreconditioning >
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static std::pair< int, double >
        call ( const OperatorImp &op,
               const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
               double reduction, double absLimit, int maxIter, bool verbose )
        {
          return solve( op.systemMatrix(), arg, dest, reduction, absLimit, maxIter, verbose );
        }

        template< class MatrixObjType, class DiscreteFunctionImp >
        static std::pair< int, double >
        solve ( const MatrixObjType &mObj,
                const DiscreteFunctionImp &arg,
                DiscreteFunctionImp &dest,
                double reduction, double absLimit, int maxIter, bool verbose )
        {
          typedef typename MatrixObjType :: MatrixAdapterType MatrixAdapterType;
          MatrixAdapterType matrix = mObj.matrixAdapter();
          
          typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;

          // verbose only in verbose mode and for rank 0 
          int verb = (verbose && (dest.space().gridPart().comm().rank() == 0)) ? 2 : 0;
            
          if( absLimit < std::numeric_limits< double >::max() )
          {
            const double residuum = matrix.residuum( arg.blockVector(), dest.blockVector() );
            reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

            if( verbose ) 
              std::cout << "ISTL CG-Solver: reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit << std::endl;
          }

          CGSolver< BlockVectorType >
            solver( matrix, matrix.scp(), matrix.preconditionAdapter(), reduction, maxIter, verb );

          // call solver and return info
          int iter = ISTLSolverCall::solve( solver, arg.blockVector(), dest.blockVector() );
          return std::pair< int, double > ( iter, matrix.averageCommTime() );
        }
      };

    public:
      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       *  \param[in] verbose   verbosity
       */
      ISTLCGOp ( const OperatorType &op,
                 double  reduction,
                 double absLimit,
                 int maxIter,
                 bool verbose )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      ISTLCGOp ( const OperatorType &op,
                 double reduction,
                 double absLimit,
                 int maxIter = std::numeric_limits< int >::max() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: ISTL CG solver,  eps = " << absLimit_ ;
        out  << "\\\\ \n";
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }


      /** \brief solve the system 
          \param[in] arg right hand side 
          \param[out] dest solution 
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        typedef SolverCaller< OperatorType, true > Caller;

        std::pair< int, double > info
          = Caller::call( op_, arg, dest, reduction_, absLimit_, maxIter_, verbose_ );

        iterations_ = info.first; 
        averageCommTime_ = info.second;
      }

      //! return number of iterations 
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
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

    private:
      // no const reference, we make const later 
      const OperatorType &op_;
      const double reduction_;
      const double absLimit_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;
      mutable double averageCommTime_;
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
    : public Operator< DF, DF >
    {
      typedef DF DiscreteFunctionType;
      typedef Op OperatorType;

    private:
      template <class OperatorImp, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static std::pair< int, double >
        call ( const OperatorImp &op,
               const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
               bool verbose )
        {
          return solve( op.systemMatrix(), arg, dest, verbose );
        }

        template< class MatrixObjType, class DiscreteFunctionImp >
        static std::pair< int, double >
        solve ( const MatrixObjType &mObj,
                const DiscreteFunctionImp &arg, DiscreteFunctionImp &dest,
                bool verbose )
        {
          // we need the type of the BCRSMatrix for the SuperLU solver
          typedef typename MatrixObjType::MatrixAdapterType MatrixAdapterType;
          typedef typename MatrixAdapterType :: MatrixType ImprovedMatrixType ;

          // we need the BCRSMatrix type here, since SuperLU is only working 
          // with that explicit type 
          typedef typename ImprovedMatrixType::BaseType MatrixType ;
          MatrixAdapterType matrix = mObj.matrixAdapter();
          
          typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;

          InverseOperatorResult returnInfo;
#if HAVE_SUPERLU
          // create solver 
          SuperLU< MatrixType > solver( matrix.getmat(), verbose );
          // solve the system 
          solver.apply( dest.blockVector(), arg.blockVector(), returnInfo );
#else 
          DUNE_THROW(NotImplemented,"SuperLU solver not found in configure, please re-configure");
#endif

          // get information 
          std::pair< int, double > p( returnInfo.iterations, matrix.averageCommTime() );
          return p; 
        }
      };

    public:
      /** \brief constructor
       *
       *  \param[in] op Mapping describing operator to invert
       *  \param[in] reduction reduction epsilon
       *  \param[in] absLimit absolut limit of residual (not used here)
       *  \param[in] maxIter maximal iteration steps
       *  \param[in] verbose verbosity
       *
       *  \note ISTL SuperLU is only available when SuperLU package is found.
       */
      ISTLSuperLU ( const OperatorType &op,
                    double  reduction,
                    double absLimit,
                    int maxIter,
                    bool verbose )
      : op_( op ),
        verbose_( verbose ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       *
       *  \note ISTL SuperLU is only available when SuperLU package is found.
       */
      ISTLSuperLU ( const OperatorType &op,
                    double reduction,
                    double absLimit,
                    int maxIter = std::numeric_limits< int >::max() )
      : op_( op ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        averageCommTime_( 0.0 )
      {}

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {}

      void finalize () const
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: ISTL SuperLU  ";
        out  << "\\\\ \n";
      }

      /** \brief solve the system 
          \param[in] arg right hand side 
          \param[out] dest solution 
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        typedef SolverCaller< OperatorType, true > Caller;

        std::pair< int, double > info
          = Caller::call( op_, arg, dest, verbose_ );

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
        apply(arg,dest);
      }

    private:
      const OperatorType &op_;
      bool verbose_ ;
      mutable int iterations_;
      mutable double averageCommTime_;
    }; 


  ///@}

  } // namespace Fem 

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: ISTLBICGSTABOp ;
  using Fem :: ISTLCGOp ;
  using Fem :: ISTLGMResOp ;
  using Fem :: ISTLSuperLU ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune 

#endif // #if HAVE_DUNE_ISTL 

#endif // #ifndef DUNE_FEM_ISTLSOLVERS_HH 
