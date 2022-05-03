#ifndef DUNE_FEM_CGINVERSEOPERATOR_HH
#define DUNE_FEM_CGINVERSEOPERATOR_HH

#include <limits>
#include <memory>
#include <type_traits>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>
#include <dune/fem/solver/parameter.hh>

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
      static_assert( (std::is_same< DomainFunctionType, RangeFunctionType >::value),
                          "DomainFunctionType must equal RangeFunctionType." );

    public:
      /** \brief constructor
       *
       *  \param[in]  epsilon        tolerance
       *  \param[in]  maxIterations  maximum number of CG iterations
       *  \param[in]  errorMeasure   use absolute (0) or relative (1) error
       *  \param[in]  verbose        verbose output
       */
      ConjugateGradientSolver ( const RealType &epsilon,
                                unsigned int maxIterations,
                                int errorMeasure,
                                bool verbose )
      : epsilon_( epsilon ),
        maxIterations_( maxIterations ),
        errorMeasure_( errorMeasure ),
        verbose_( verbose && Fem::Parameter::verbose() ),
        averageCommTime_( 0.0 ),
        realCount_( 0 )
      {}
      ConjugateGradientSolver ( const RealType &epsilon,
                                unsigned int maxIterations,
                                bool verbose,
                                const ParameterReader &parameter = Parameter::container() )
      : epsilon_( epsilon ),
        maxIterations_( maxIterations ),
        errorMeasure_( 1 ),
        verbose_( verbose && Fem::Parameter::verbose() ),
        averageCommTime_( 0.0 ),
        realCount_( 0 )
      {}

      /** \brief constructor
       *
       *  \param[in]  epsilon        tolerance
       *  \param[in]  maxIterations  maximum number of CG iterations
       */
      ConjugateGradientSolver ( RealType epsilon,
                                unsigned int maxIterations,
                                const ParameterReader &parameter = Parameter::container() )
      : epsilon_( epsilon ),
        maxIterations_( maxIterations ),
        errorMeasure_( 1 ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        averageCommTime_( 0.0 ),
        realCount_( 0 )
      {}

      // ConjugateGradientSolver ( const ThisType & )=delete;

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

      void setMaxIterations ( unsigned int maxIterations ) { maxIterations_ = maxIterations; }


      //! return average communication time during last solve
      double averageCommTime() const
      {
        return averageCommTime_;
      }

    protected:
      const RealType epsilon_;
      unsigned int maxIterations_;
      int errorMeasure_;
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
        typedef CGInverseOperator< DiscreteFunction >               ThisType;

      public:
        typedef typename BaseType::DomainFunctionType DomainFunctionType;
        typedef typename BaseType::RangeFunctionType RangeFunctionType;

        typedef Fem::Operator< DomainFunctionType, RangeFunctionType > OperatorType;
        typedef Fem::Operator< RangeFunctionType, DomainFunctionType > PreconditionerType;

        typedef typename OperatorType::RangeFieldType RangeFieldType;
        typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

        // non-const version
        using BaseType::finalize;

      private:
        typedef ConjugateGradientSolver< OperatorType > SolverType;

      public:
        /** \brief constructor of CGInverseOperator
         *
         *  \param[in]  redEps   reduction epsilon
         *  \param[in]  absLimit absolut limit of residual
         *  \param[in]  maxIter  maximum number of iteration steps
         *  \param[in]  verbose  verbosity
         */
        CGInverseOperator ( RealType redEps, RealType absLimit,
                            unsigned int maxIter, bool verbose,
                            const ParameterReader &parameter = Parameter::container() )
          : solver_( absLimit, maxIter, verbose, parameter ),
            parameter_( parameter )
        {}

        CGInverseOperator ( const SolverParameter& param = SolverParameter(Parameter::container() ) )
          : solver_( param.tolerance(),
                     param.maxIterations(),
                     param.errorMeasure(),
                     param.verbose() ),
            parameter_( param )
        {
        }


        /** \brief constructor of CGInverseOperator
         *
         *  \param[in]  redEps   reduction epsilon
         *  \param[in]  absLimit absolut limit of residual
         *  \param[in]  maxIter  maximum number of iteration steps
         */
        CGInverseOperator ( RealType redEps, RealType absLimit,
                            unsigned int maxIter,
                            const ParameterReader &parameter = Parameter::container() )
          : solver_( absLimit, maxIter, parameter ),
            parameter_( parameter )
        {}

        CGInverseOperator ( RealType redEps, RealType absLimit,
                            const ParameterReader &parameter = Parameter::container() )
          : CGInverseOperator( redEps, absLimit, std::numeric_limits< unsigned int >::max(), parameter )
        {}

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
                            unsigned int maxIter, bool verbose,
                            const ParameterReader &parameter = Parameter::container() )
          : CGInverseOperator( redEps, absLimit, maxIter, verbose, parameter )
        {
          bind( op );
        }


        /** \brief constructor of CGInverseOperator
         *
         *  \param[in]  op       operator to invert
         *  \param[in]  redEps   reduction epsilon
         *  \param[in]  absLimit absolut limit of residual
         *  \param[in]  maxIter  maximum number of iteration steps
         */
        CGInverseOperator ( const OperatorType &op,
                            RealType redEps, RealType absLimit,
                            unsigned int maxIter,
                            const ParameterReader &parameter = Parameter::container() )
          : CGInverseOperator( redEps, absLimit, maxIter, parameter )
        {
          bind( op );
        }

        CGInverseOperator ( const OperatorType &op,
                            RealType redEps, RealType absLimit,
                            const ParameterReader &parameter = Parameter::container() )
          : CGInverseOperator( redEps, absLimit, parameter )
        {
          bind( op );
        }

        /** \brief constructor of CGInverseOperator
         *
         *  \param[in]  op       operator to invert
         *  \param[in]  precond  precondition operator
         *  \param[in]  redEps   reduction epsilon
         *  \param[in]  absLimit absolut limit of residual
         *  \param[in]  maxIter  maximum number of iteration steps
         *  \param[in]  verbose  verbosity
         */
        CGInverseOperator ( const OperatorType &op,
                            const PreconditionerType &precond,
                            RealType redEps, RealType absLimit,
                            unsigned int maxIter, bool verbose,
                            const ParameterReader &parameter = Parameter::container() )
          : CGInverseOperator( redEps, absLimit, maxIter, verbose )
        {
          bind( op, precond );
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
                            const PreconditionerType &precond,
                            RealType redEps, RealType absLimit,
                            const ParameterReader &parameter = Parameter::container() )
          : CGInverseOperator( redEps, absLimit, parameter )
        {
          bind( op, precond );
        }

        CGInverseOperator ( const OperatorType &op,
                            const PreconditionerType &precond,
                            RealType redEps, RealType absLimit,
                            unsigned int maxIter,
                            const ParameterReader &parameter = Parameter::container() )
          : CGInverseOperator( redEps, absLimit, maxIter, parameter )
        {
          bind( op, precond );
        }

        void bind ( const OperatorType &op ) { operator_ = &op; preconditioner_ = nullptr; }
        void bind ( const OperatorType &op, const PreconditionerType &precond )
        {
          operator_ = &op;
          preconditioner_ = &precond;
        }
        void unbind () { operator_ = nullptr; preconditioner_ = nullptr; }

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
          const_cast< ThisType& > (*this).finalize();
        }

        template<typename... A>
        inline void prepare(A... ) const
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
          assert(operator_);
          if(preconditioner_)
            solver_.solve( *operator_, *preconditioner_, arg, dest );
          else
            solver_.solve( *operator_, arg, dest );
        }

        //! number of iterations needed for last solve
        unsigned int iterations () const
        {
          return solver_.iterations();
        }

        void setMaxIterations ( unsigned int maxIter ) { solver_.setMaxIterations( maxIter ); }

        //! return average communication time during last solve
        double averageCommTime() const
        {
          return solver_.averageCommTime();
        }

      SolverParameter& parameter () const
      {
        return parameter_;
      }

      protected:
        const OperatorType *operator_ = nullptr;
        const PreconditionerType *preconditioner_ = nullptr;
        SolverType solver_;
        mutable SolverParameter parameter_;
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
      using BaseType::bind;

      typedef SolverParameter SolverParameterType;
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef DomainFunctionType DestinationType;

      //! type of operator
      typedef Op OperatorType;

      typedef typename OperatorType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

      // Preconditioner is to approximate op^-1 !
      typedef Fem::Operator< RangeFunctionType,DomainFunctionType > PreconditioningType;

      CGInverseOperator ( const SolverParameter& param = SolverParameter(Parameter::container()) )
        : BaseType( param )
      {}


      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       *  \param[in]  verbose  verbosity
       */
      CGInverseOperator ( RealType redEps, RealType absLimit,
                          unsigned int maxIter, bool verbose,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( redEps, absLimit, maxIter, verbose, parameter )
      {}

      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       */
      CGInverseOperator ( RealType redEps, RealType absLimit,
                          unsigned int maxIter,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( redEps, absLimit, maxIter, parameter )
      {}

      CGInverseOperator ( RealType redEps, RealType absLimit,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( redEps, absLimit, parameter )
      {}

      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       *  \param[in]  verbose  verbosity
       */
      template< class LinearOperator, std::enable_if_t< std::is_base_of< OperatorType, LinearOperator >::value, int > = 0 >
      CGInverseOperator ( const LinearOperator &op,
                          RealType redEps, RealType absLimit,
                          unsigned int maxIter, bool verbose,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( redEps, absLimit, maxIter, verbose, parameter )
      {
        bind( op );
      }

      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       */
      template< class LinearOperator, std::enable_if_t< std::is_base_of< OperatorType, LinearOperator >::value, int > = 0 >
      CGInverseOperator ( const LinearOperator &op,
                          RealType redEps, RealType absLimit,
                          unsigned int maxIter,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( redEps, absLimit, maxIter, parameter )
      {
        bind( op );
      }

      template< class LinearOperator, std::enable_if_t< std::is_base_of< OperatorType, LinearOperator >::value, int > = 0 >
      CGInverseOperator ( const LinearOperator &op,
                          RealType redEps, RealType absLimit,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( redEps, absLimit, parameter )
      {
        bind( op );
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
                          unsigned int maxIter, bool verbose,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( op, precond, redEps, absLimit, maxIter, verbose, parameter )
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
                          const PreconditioningType &precond,
                          RealType redEps, RealType absLimit,
                          unsigned int maxIter,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( op, precond, redEps, absLimit, maxIter, parameter )
      {}

      CGInverseOperator ( const OperatorType &op,
                          const PreconditioningType &precond,
                          RealType redEps, RealType absLimit,
                          const ParameterReader &parameter = Parameter::container() )
      : BaseType( op, precond, redEps, absLimit, parameter )
      {}


      template< class LinearOperator, std::enable_if_t< std::is_base_of< OperatorType, LinearOperator >::value, int > = 0 >
      void bind ( const LinearOperator &op )
      {
        BaseType::bind( op );
        checkPreconditioning( op );
      }

      void unbind ()
      {
        BaseType::unbind();
        precondObj_.reset();
      }

    protected:
      template< class LinearOperator >
      void checkPreconditioning( const LinearOperator &linearOp )
      {
        bool preconditioning = false;
        if (!parameter_.parameter().exists(parameter_.keyPrefix()+"preconditioning.method") &&
            parameter_.parameter().exists("fem.preconditioning"))
        {
          preconditioning = parameter_.parameter().template getValue< bool >( "fem.preconditioning", false );
          std::cout << "WARNING: using deprecated parameter `fem.preconditioning` use "
                    << parameter_.keyPrefix() << "preconditioning.method instead\n";
        }
        else
          preconditioning = parameter_.preconditionMethod(
                {
                  SolverParameter::none,   // no preconditioning
                  SolverParameter::jacobi  // Jacobi preconditioning
                });
        if( preconditioning && std::is_base_of< AssembledOperator< DomainFunctionType, DomainFunctionType > ,LinearOperator > :: value )
        {
          // create diagonal preconditioner
          precondObj_.reset( new DiagonalPreconditioner< DomainFunctionType, LinearOperator >( linearOp ) );
          preconditioner_ = precondObj_.get();
        }
      }

      using BaseType::preconditioner_;
      using BaseType::parameter_;
      std::unique_ptr< PreconditioningType > precondObj_;
    };

    // Implementation of ConjugateGradientSolver
    // -----------------------------------------

    template< class Operator >
    inline void ConjugateGradientSolver< Operator >
      ::solve ( const OperatorType &op, const RangeFunctionType &b, DomainFunctionType &x ) const
    {

      RealType tolerance = (epsilon_ * epsilon_);
      if (errorMeasure_ == 1)
        tolerance *= b.normSquaredDofs( );

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
        if( verbose_ )
        {
          std::cout << "CG-Iteration: " << realCount_ << ", Residuum: " << std::sqrt(residuum)
            << ", Tolerance: " << std::sqrt(tolerance) << std::endl;
          // only for parallel apps
          if( b.space().gridPart().comm().size() > 1 )
            std::cout << "Communication needed: " << exchangeTime << " s" << std::endl;
        }

        averageCommTime_ += exchangeTime;
      }
    }


    template< class Operator >
    inline void ConjugateGradientSolver< Operator >
    ::solve ( const OperatorType &op, const PreconditionerType &precond, const RangeFunctionType &b, DomainFunctionType &x ) const
    {
      RealType tolerance = (epsilon_ * epsilon_);
      if (errorMeasure_ == 1)
        tolerance *= b.normSquaredDofs( );

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
        if( verbose_ )
        {
          std::cout << "CG-Iteration: " << realCount_ << ", Residuum: " << std::sqrt(residuum)
            << ", Tolerance: " << std::sqrt(tolerance) << std::endl;
          // only for parallel apps
          if( b.space().gridPart().comm().size() > 1 )
            std::cout << "Communication needed: " << exchangeTime << " s" << std::endl;
        }

        averageCommTime_ += exchangeTime;
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_CGINVERSEOPERATOR_HH
