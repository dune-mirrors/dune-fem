#ifndef DUNE_FEM_INVERSEOPERATORS_HH
#define DUNE_FEM_INVERSEOPERATORS_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  // ConjugateGradientSolver
  // -----------------------

  /** \class ConjugateGradientSolver
   *  \ingroup OEMSolver
   *  \brief   linear solver using the CG algorithm
   *
   *  \param  Operator  type of the operator to invert
   */
  template< class Operator >
  class ConjugateGradientSolver
  {
    typedef ConjugateGradientSolver< Operator > ThisType;

  public:
    //! type of the operators to invert
    typedef Operator OperatorType;

    //! field type of the operator's domain vectors
    typedef typename OperatorType::DomainFieldType DomainFieldType;
    //! field type of the operator's range vectors
    typedef typename OperatorType::RangeFieldType RangeFieldType;
    
    //! type of the operator's domain vectors
    typedef typename OperatorType::DomainFunctionType DomainFunctionType;
    //! type of the operator's range vectors
    typedef typename OperatorType::RangeFunctionType RangeFunctionType;

  private:
    dune_static_assert( (Conversion< DomainFunctionType, RangeFunctionType >::sameType),
                        "DomainFunctionType must equal RangeFunctionType." );

  public:
    /** \brief constructor
     *
     *  \param[in]  epsilon        tolerance
     *  \param[in]  maxIterations  maximum number of CG iterations
     *  \param[in]  verbose        verbose output
     */
    ConjugateGradientSolver ( const RangeFieldType &epsilon,
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
    ConjugateGradientSolver ( RangeFieldType epsilon,
                              unsigned int maxIterations )
    : epsilon_( epsilon ),
      maxIterations_( maxIterations ),
      verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
      averageCommTime_( 0.0 ),
      realCount_( 0 )
    {}

  private:
    // prohibit copying
    ConjugateGradientSolver ( const ThisType & );

  public:
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
    const RangeFieldType epsilon_;
    const unsigned int maxIterations_;
    const bool verbose_;
    mutable double averageCommTime_;
    mutable unsigned int realCount_;
  };



  // CGInverseOperator
  // -----------------
 
  /** \class   CGInverseOperator
   *  \ingroup OEMSolver
   *  \brief   Inverse operator base on CG method
   */
  template< class DiscreteFunction >
  class CGInverseOperator
  : public Fem::Operator< DiscreteFunction, DiscreteFunction >
  {
    typedef Fem::Operator< DiscreteFunction, DiscreteFunction > BaseType;
    
  public:
    typedef typename BaseType::DomainFunctionType DomainFunctionType;
    typedef typename BaseType::RangeFunctionType RangeFunctionType;

    typedef Fem::Operator< DiscreteFunction, DiscreteFunction > OperatorType;

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
                        double redEps, double absLimit,
                        unsigned int maxIter, bool verbose )
    : operator_( op ),
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
                        double redEps, double absLimit,
                        unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
    : operator_( op ),
      solver_( absLimit, maxIter )
    {}

    /** \brief application operator
     *
     *  The application operator actually solves the linear system
     *  \f$op(w) = u\f$ using the CG method.
     *
     *  \param[in]   u  argument discrete function
     *  \param[out]  w  destination discrete function
     */
    virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      solver_.solve( operator_, u, w );
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
    SolverType solver_;
  };



  /** \class   CGInverseOp
   *  \ingroup OEMSolver
   *  \brief   Inversion operator using CG algorithm, operator type is a
   *           template parameter.
   */
  template< class DF, class Op >
  struct CGInverseOp
  : public Operator< typename DF::RangeFieldType, typename DF::RangeFieldType, DF, DF >
  {
    typedef DF DiscreteFunctionType;
    typedef Op OperatorType;

    /** \brief constructor of CGInverseOperator
     *
     *  \param[in] op Mapping describing operator to invert
     *  \param[in] redEps reduction epsilon
     *  \param[in] absLimit absolut limit of residual
     *  \param[in] maxIter maximal iteration steps
     *  \param[in] verbose verbosity
     */
    CGInverseOp( const OperatorType &op,
                 double redEps,
                 double absLimit,
                 int maxIter,
                 bool verbose )
    : operator_( op ),
      solver_( absLimit, maxIter, verbose )
    {} 

    /** \brief constructor of CGInverseOperator
     *
     *  \param[in] op Mapping describing operator to invert
     *  \param[in] redEps reduction epsilon
     *  \param[in] absLimit absolut limit of residual
     *  \param[in] maxIter maximal iteration steps
     */
    CGInverseOp( const OperatorType &op,
                 double redEps,
                 double absLimit,
                 unsigned int maxIter = std::numeric_limits< int >::max() )
    : operator_( op ),
      solver_( absLimit, maxIter )
    {} 

    /** \brief solve the system 
        \param[in] arg right hand side 
        \param[out] dest solution 
    */
    virtual void operator() ( const DiscreteFunctionType &arg,
                              DiscreteFunctionType &dest ) const
    {
      solver_.solve( operator_, arg, dest );
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
    const ConjugateGradientSolver< OperatorType > solver_;
  };



  // Implementation of ConjugateGradientSolver
  // -----------------------------------------

  template< class Operator >
  inline void ConjugateGradientSolver< Operator >
    ::solve ( const OperatorType &op, const RangeFunctionType &b, DomainFunctionType &x ) const
  {
    const bool verbose = (verbose_ && (b.space().grid().comm().rank() == 0));
    
    const RangeFieldType tolerance = (epsilon_ * epsilon_) * b.scalarProductDofs( b ); 

    averageCommTime_ = 0.0;
    
    RangeFunctionType h( b );
    op( x, h );

    RangeFunctionType r( h );
    r -= b;

    RangeFunctionType p( b );
    p -= h;

    RangeFieldType prevResiduum = 0;
    RangeFieldType residuum = r.scalarProductDofs( r );
 
    for( realCount_ = 0; (residuum > tolerance) && (realCount_ < maxIterations_); ++realCount_ )
    {
      if( realCount_ > 0 )
      { 
        p *= (residuum / prevResiduum);
        p -= r;
      }

      op( p, h );

      const RangeFieldType alpha = residuum / p.scalarProductDofs( h );
      x.addScaled( p, alpha );
      r.addScaled( h, alpha );

      prevResiduum = residuum;
      residuum = r.scalarProductDofs( r );
      
      double exchangeTime = h.space().communicator().exchangeTime();
      if( verbose )
      {
        std::cerr << "CG-Iteration: " << realCount_ << ", Residuum: " << residuum << std::endl;
        // only for parallel apps 
        if( b.space().grid().comm().size() > 1 )
          std::cerr << "Communication needed: " << exchangeTime << " s" << std::endl;
      }
      
      averageCommTime_ += exchangeTime;
    }
  }

} // end namespace Dune

#endif // #ifndef DUNE_FEM_INVERSEOPERATORS_HH
