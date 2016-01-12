#ifndef DUNE_FEM_CGINVERSEOPERATOR_HH
#define DUNE_FEM_CGINVERSEOPERATOR_HH

#include <type_traits>

#include <dune/common/typetraits.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>

#include <dune/fem/solver/diagonalpreconditioner.hh>

namespace Dune
{

  namespace Fem
  {

    // ConjugateGradientSolver
    // -----------------------

    /** \class ConjugateGradientSolver
     *  \ingroup OEMSolver
     *  \brief   linear solver using the CG algorithm
     *
     *  \param  Operator  type of the operator to invert
     */
    template< class Operator>
    class ConjugateGradientSolver
    {
      typedef ConjugateGradientSolver< Operator> ThisType;

    public:
      //! type of the operators to invert
      typedef Operator OperatorType;

      //! field type of the operator's domain vectors
      typedef typename OperatorType::DomainFieldType DomainFieldType;
      //! field type of the operator's range vectors
      typedef typename OperatorType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

      //! type of the operator's domain vectors
      typedef typename OperatorType::DomainFunctionType DomainFunctionType;
      //! type of the operator's range vectors
      typedef typename OperatorType::RangeFunctionType RangeFunctionType;

      //! type of the preconditioner, maps from the range of the operator (the dual space) in it's domain
      typedef Fem::Operator< RangeFunctionType, DomainFunctionType > PreconditionerType;


    private:
      static_assert( (Conversion< DomainFunctionType, RangeFunctionType >::sameType),
                          "DomainFunctionType must equal RangeFunctionType." );

    public:
      /** \brief constructor
       *
       *  \param[in]  epsilon        tolerance
       *  \param[in]  maxIterations  maximum number of CG iterations
       *  \param[in]  verbose        verbose output
       */
      ConjugateGradientSolver ( const RealType &epsilon,
                                unsigned int maxIterations,
                                bool verbose )
      : epsilon_( epsilon ),
        maxIterations_( maxIterations ),
        verbose_( verbose ),
        averageCommTime_( 0.0 ),
        realCount_( 0 )

      {}

      /** \brief constructor
       *
       *  \param[in]  epsilon        tolerance
       *  \param[in]  maxIterations  maximum number of CG iterations
       */
      ConjugateGradientSolver ( RealType epsilon,
                                unsigned int maxIterations )
      : epsilon_( epsilon ),
        maxIterations_( maxIterations ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        averageCommTime_( 0.0 ),
        realCount_( 0 )
      {}

      ConjugateGradientSolver ( const ThisType & )=delete;

      /** \brief solve \f$op( x ) = b\f$
       *
       *  \note The CG algorithm also works for positive semidefinite operators.
       *        In this case, \f$x \cdot v = b \cdot v\f$ for all \f$v\f$ in the
       *        operator's kernel.
       *
       *  \param[in]   op  linear operator to invert (must be symmetic and
       *                   positive definite)
       *  \param[in]   b   right hand side
       *  \param       x   solution (must be initialized to a start value)
       */
      void solve ( const OperatorType &op, const RangeFunctionType &b, DomainFunctionType &x ) const;

      /** \brief solve \f$op( x ) = b\f$
       *
       *  \note The CG algorithm also works for positive semidefinite operators.
       *        In this case, \f$x \cdot v = b \cdot v\f$ for all \f$v\f$ in the
       *        operator's kernel.
       *
       *  \param[in]   op  linear operator to invert (must be symmetic and
       *                   positive definite)
       *  \param[in]   p   (lef) preconditioning operator
       *  \param[in]   b   right hand side
       *  \param       x   solution (must be initialized to a start value)
       */
      void solve ( const OperatorType &op, const PreconditionerType &p,
                   const RangeFunctionType &b, DomainFunctionType &x ) const;

      //! number of iterations needed for last solve
      unsigned int iterations () const
      {
        return realCount_;
      }

      //! return average communication time during last solve
      double averageCommTime() const
      {
        return averageCommTime_;
      }

    protected:
      const RealType epsilon_;
      const unsigned int maxIterations_;
      const bool verbose_;
      mutable double averageCommTime_;
      mutable unsigned int realCount_;
    };


    namespace Solver
    {
      // CGInverseOperator
      // -----------------

      /** \class   CGInverseOperator
       *  \ingroup OEMSolver
       *  \brief   Inverse operator base on CG method. This is the base class for the
       *           cg solver and does not imvolve any runtime parametrization
       */
      template< class DiscreteFunction >
      class CGInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
      {
        typedef Fem::Operator< DiscreteFunction, DiscreteFunction > BaseType;

      public:
        typedef typename BaseType::DomainFunctionType DomainFunctionType;
        typedef typename BaseType::RangeFunctionType RangeFunctionType;

        typedef Fem::Operator< DomainFunctionType, RangeFunctionType > OperatorType;
        typedef Fem::Operator< RangeFunctionType, DomainFunctionType > PreconditionerType;

        typedef typename OperatorType::RangeFieldType RangeFieldType;
        typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

      private:
        typedef ConjugateGradientSolver< OperatorType > SolverType;

      public:
        /** \brief constructor of CGInverseOperator
         *
         *  \param[in]  op       operator to invert
         *  \param[in]  redEps   reduction epsilon
         *  \param[in]  absLimit absolut limit of residual
         *  \param[in]  maxIter  maximum number of iteration steps
         *  \param[in]  verbose  verbosity
         */
        CGInverseOperator ( const OperatorType &op,
                            RealType redEps, RealType absLimit,
                            unsigned int maxIter, bool verbose )
          : operator_( op ),
            preconditioner_ ( 0 ),
            solver_( absLimit, maxIter, verbose )
        {}

        /** \brief constructor of CGInverseOperator
         *
         *  \param[in]  op       operator to invert
         *  \param[in]  redEps   reduction epsilon
         *  \param[in]  absLimit absolut limit of residual
         *  \param[in]  maxIter  maximum number of iteration steps
         */
        CGInverseOperator ( const OperatorType &op,
                            RealType redEps, RealType absLimit,
                            unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
          : operator_( op ),
            preconditioner_ ( 0 ),
            solver_( absLimit, maxIter )
        {}

        /** \brief constructor of CGInverseOperator
         *
         *  \param[in]  op       operator to invert
         *  \param[in]  precond  precondition operator
         *  \param[in]  redEps   reduction epsilon
         *  \param[in]  absLimit absolut limit of residual
         *  \param[in]  maxIter  maximum number of iteration steps
         */
        CGInverseOperator ( const OperatorType &op,
                            const PreconditionerType &precond,
                            RealType redEps, RealType absLimit,
                            unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
          : operator_( op ),
            preconditioner_( &precond ),
            solver_( absLimit, maxIter )
        {}

        /** \brief application operator
         *
         *  The application operator actually solves the linear system
         *  \f$op(dest) = arg\f$ using the CG method.
         *
         *  \param[in]   arg  argument discrete function
         *  \param[out]  dest  destination discrete function
         */
        virtual void operator()( const DomainFunctionType &arg, RangeFunctionType &dest ) const
        {
          prepare();
          apply(arg,dest);
          finalize();
        }

        template<typename... A>
        inline void prepare(A... ) const
        {}

        inline void finalize() const
        {}

        /** \brief application operator
         *
         *  The application operator actually solves the linear system
         *  \f$op(dest) = arg\f$ using the CG method.
         *
         *  \param[in]   arg  argument discrete function
         *  \param[out]  dest  destination discrete function
         */
        virtual void apply( const DomainFunctionType &arg, RangeFunctionType &dest ) const
        {
          if(preconditioner_)
            solver_.solve( operator_, *preconditioner_, arg, dest );
          else
            solver_.solve(operator_,arg,dest);
        }

        //! number of iterations needed for last solve
        unsigned int iterations () const
        {
          return solver_.iterations();
        }

        //! return average communication time during last solve
        double averageCommTime() const
        {
          return solver_.averageCommTime();
        }

      protected:
        const OperatorType &operator_;
        const PreconditionerType *preconditioner_;
        SolverType solver_;
      };
    }

    // CGInverseOperator
    // -----------------

    /** \class   CGInverseOperator
     *  \ingroup OEMSolver
     *  \brief   Inverse operator base on CG method. Uses a runtime parameter
     *           fem.preconditioning which enables diagonal preconditioning if
     *           diagonal matrix entries are available, i.e.,
     *           Op :: assembled is true.
     */
    template< class DiscreteFunction,
              class Op = Fem::Operator< DiscreteFunction, DiscreteFunction > >
    class CGInverseOperator
    : public Fem::Solver::CGInverseOperator< DiscreteFunction >
    {
      typedef Fem::Solver::CGInverseOperator< DiscreteFunction > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef DomainFunctionType DestinationType;

      //! type of operator
      typedef Op OperatorType;

      typedef typename OperatorType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

      // Preconditioner is to approximate op^-1 !
      typedef Fem::Operator< RangeFunctionType,DomainFunctionType > PreconditioningType;

      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       *  \param[in]  verbose  verbosity
       */
      template <class LinearOperator>
      CGInverseOperator ( const LinearOperator &op,
                          RealType redEps, RealType absLimit,
                          unsigned int maxIter, bool verbose )
      : BaseType( op, redEps, absLimit, maxIter, verbose ),
        precondObj_( 0 )
      {
        checkPreconditioning( op );
      }

      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       */
      template <class LinearOperator>
      CGInverseOperator ( const LinearOperator &op,
                          RealType redEps, RealType absLimit,
                          unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
      : BaseType( op, redEps, absLimit, maxIter ),
        precondObj_( 0 )
      {
        checkPreconditioning( op );
      }

      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  precond  precondition operator
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       */
      CGInverseOperator ( const OperatorType &op,
                          const PreconditioningType &precond,
                          RealType redEps, RealType absLimit,
                          unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
      : BaseType( op, precond, redEps, absLimit, maxIter ),
        precondObj_( 0 )
      {}

      /** \brief destructor */
      ~CGInverseOperator ()
      {
        if( precondObj_ )
          delete precondObj_;
      }

    protected:
      template< class LinearOperator >
      void checkPreconditioning( const LinearOperator &linearOp )
      {
        const bool preconditioning = Parameter::getValue< bool >( "fem.preconditioning", false );
        if( preconditioning && std::is_base_of< AssembledOperator< DomainFunctionType, DomainFunctionType > ,LinearOperator > :: value )
        {
          // create diagonal preconditioner
          precondObj_ = new DiagonalPreconditioner< DomainFunctionType, LinearOperator >( linearOp );
          preconditioner_ = precondObj_;
        }
      }

      using BaseType::preconditioner_;
      PreconditioningType *precondObj_;
    };

    // Implementation of ConjugateGradientSolver
    // -----------------------------------------

    template< class Operator >
    inline void ConjugateGradientSolver< Operator >
      ::solve ( const OperatorType &op, const RangeFunctionType &b, DomainFunctionType &x ) const
    {
      const bool verbose = (verbose_ && (b.space().gridPart().comm().rank() == 0));

      const RealType tolerance = (epsilon_ * epsilon_) * b.normSquaredDofs( );

      averageCommTime_ = 0.0;

      RangeFunctionType h( b );
      op( x, h );

      RangeFunctionType r( h );
      r -= b;

      RangeFunctionType p( b );
      p -= h;

      RealType prevResiduum = 0;
      RealType residuum = r.normSquaredDofs( );

      for( realCount_ = 0; (residuum > tolerance) && (realCount_ < maxIterations_); ++realCount_ )
      {
        if( realCount_ > 0 )
        {
          assert( residuum/prevResiduum == residuum/prevResiduum );
          p *= (residuum / prevResiduum);
          p -= r;
        }

        op( p, h );

        RangeFieldType pdoth = p.scalarProductDofs( h );
        const RangeFieldType alpha = residuum / pdoth;
        assert( alpha == alpha );
        x.axpy( alpha, p );
        r.axpy( alpha, h );

        prevResiduum = residuum;
        residuum = r.normSquaredDofs( );

        double exchangeTime = h.space().communicator().exchangeTime();
        if( verbose )
        {
          std::cerr << "CG-Iteration: " << realCount_ << ", sqr(Residuum): " << residuum << std::endl;
          // only for parallel apps
          if( b.space().gridPart().comm().size() > 1 )
            std::cerr << "Communication needed: " << exchangeTime << " s" << std::endl;
        }

        averageCommTime_ += exchangeTime;
      }
    }


    template< class Operator >
    inline void ConjugateGradientSolver< Operator >
    ::solve ( const OperatorType &op, const PreconditionerType &precond, const RangeFunctionType &b, DomainFunctionType &x ) const
    {
      const bool verbose = (verbose_ && (b.space().gridPart().comm().rank() == 0));

      const RealType tolerance = (epsilon_ * epsilon_) * b.normSquaredDofs( );

      averageCommTime_ = 0.0;

      RangeFunctionType h( b );
      //h=Ax
      op( x, h );

      //r=Ax-b
      RangeFunctionType r( h );
      r -= b;

      //p=b-A*x <= r_0 Deufelhard
      RangeFunctionType p( b );
      p -= h;

      //q=B*p <=q Deuf
      RangeFunctionType q ( b );
      precond(p,q);

      RangeFunctionType s (q);

      RangeFieldType prevResiduum = 0;    // note that these will be real_type but require scalar product evaluation
      RangeFieldType residuum = p.scalarProductDofs( q );//<p,Bp>

      for( realCount_ = 0; (std::real(residuum) > tolerance) && (realCount_ < maxIterations_); ++realCount_ )
      {
        if( realCount_ > 0 )
        {
          assert( residuum/prevResiduum == residuum/prevResiduum );
          const RangeFieldType beta=residuum/prevResiduum;
          q*=beta;
          q+=(s);
        }

        op( q, h );

        RangeFieldType qdoth = q.scalarProductDofs( h );
        const RangeFieldType alpha = residuum / qdoth;//<p,Bp>/<q,Aq>
        assert( alpha == alpha );
        x.axpy( alpha, q );

        p.axpy( -alpha, h );//r_k

        precond(p,s); //B*r_k

        prevResiduum = residuum;//<rk-1,B*rk-1>

        residuum = p.scalarProductDofs( s );//<rk,B*rk>

        double exchangeTime = h.space().communicator().exchangeTime();
        if( verbose )
        {
          std::cerr << "CG-Iteration: " << realCount_ << ", Residuum: " << residuum << std::endl;
          // only for parallel apps
          if( b.space().gridPart().comm().size() > 1 )
            std::cerr << "Communication needed: " << exchangeTime << " s" << std::endl;
        }

        averageCommTime_ += exchangeTime;
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_CGINVERSEOPERATOR_HH
