#ifndef DUNE_FEM_INVERSEOPERATORS_HH
#define DUNE_FEM_INVERSEOPERATORS_HH

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  /** \class ConjugateGradientSolver
   *  \ingroup OEMSolver
   *  \brief   linear solver using the CG algorithm
   *
   *  \param  Operator  type of the operator to invert
   */
  template< class Operator >
  class ConjugateGradientSolver
  {
  public:
    //! type of the operators to invert
    typedef Operator OperatorType;

    //! field type of the operator's domain vectors
    typedef typename OperatorType :: DomainFieldType DomainFieldType;
    //! field type of the operator's range vectors
    typedef typename OperatorType :: RangeFieldType RangeFieldType;
    
    //! type of the operator's domain vectors
    typedef typename OperatorType :: DomainType DomainType;
    //! type of the operator's range vectors
    typedef typename OperatorType :: RangeType RangeType;

  private:
    typedef CompileTimeChecker< Conversion< DomainType, RangeType > :: sameType >
      __DOMAINTYPE_MUST_EQUAL_RANGETYPE__;

  protected:
    const RangeFieldType epsilon_;
    const unsigned int maxIterations_;
    const bool verbose_;
    
  public:
    /** \brief constructor
     *
     *  \param[in]  epsilon        tolerance
     *  \param[in]  maxIterations  maximum number of CG iterations
     *  \param[in]  verbose        verbose output
     */
    inline ConjugateGradientSolver ( RangeFieldType epsilon,
                                     unsigned int maxIterations,
                                     bool verbose = false )
    : epsilon_( epsilon ),
      maxIterations_( maxIterations ),
      verbose_( verbose )
    {}

  private:
    // prohibit copying
    ConjugateGradientSolver ( const ConjugateGradientSolver & );

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
    inline void solve ( const OperatorType &op,
                        const RangeType &b,
                        DomainType &x ) const
    {
      const bool verbose = (verbose_ && (b.space().grid().comm().rank() == 0));

      const RangeFieldType tolerance = SQR( epsilon_ ) * b.scalarProductDofs( b );

      RangeType h( b );
      op( x, h );

      RangeType r( h );
      r -= b;

      RangeType p( b );
      p -= h;

      RangeFieldType prevResiduum = 0;
      RangeFieldType residuum = r.scalarProductDofs( r );
   
      for( unsigned int count = 0;
           (residuum > tolerance) && (count < maxIterations_); ++count )
      {
        if( count > 0 )
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
        
        if( verbose )
          std :: cerr << "CG-Iteration: " << count << ", Residuum: " << residuum
                      << std :: endl;
      }
    }
  };


 
  /** \class   CGInverseOperator
   *  \ingroup OEMSolver
   *  \brief   Inversion operator using CG algorithm, operator type is Mapping
   */
  template <class DiscreteFunctionType>
  class CGInverseOperator : public Operator<
    typename DiscreteFunctionType::DomainFieldType,
    typename DiscreteFunctionType::RangeFieldType,
    DiscreteFunctionType,DiscreteFunctionType> 
  {
    typedef Operator<
            typename DiscreteFunctionType::DomainFieldType,
            typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> BaseType;
    
    typedef Mapping<typename DiscreteFunctionType::DomainFieldType ,
                    typename DiscreteFunctionType::RangeFieldType ,
                    DiscreteFunctionType,DiscreteFunctionType> MappingType;

  public:
    typedef typename BaseType :: DomainType DomainType;
    typedef typename BaseType :: RangeType RangeType;

  protected:
    const MappingType &operator_;

    const ConjugateGradientSolver< MappingType > solver_;

  public:
    /** \brief constructor of CGInverseOperator
      \param[in] op Mapping describing operator to invert 
      \param[in] redEps reduction epsilon 
      \param[in] absLimit absolut limit of residual  
      \param[in] maxIter maximal iteration steps 
      \param[in] verbose verbosity 
    */
    CGInverseOperator( const MappingType &op,
                       double redEps,
                       double absLimit,
                       int maxIter,
                       int verbose )
    : operator_( op ),
      solver_( absLimit, maxIter, (verbose > 0) )
    {}

    /** \brief solve the system 
        \param[in] arg right hand side 
        \param[out] dest solution 
    */
    virtual void operator() ( const DomainType &arg,
                              RangeType &dest ) const
    {
      solver_.solve( operator_, arg, dest );
    }
  };



  /** \class   CGInverseOp
   *  \ingroup OEMSolver
   *  \brief   Inversion operator using CG algorithm, operator type is a
   *           template parameter.
   */
  template <class DiscreteFunctionType, class OperatorType>
  class CGInverseOp : public Operator<
    typename DiscreteFunctionType::DomainFieldType,
    typename DiscreteFunctionType::RangeFieldType,
    DiscreteFunctionType,DiscreteFunctionType> 
  {
  protected:
    const OperatorType &operator_;

    const ConjugateGradientSolver< OperatorType > solver_;

  public:
    /** \brief constructor of CGInverseOperator
      \param[in] op Mapping describing operator to invert 
      \param[in] redEps reduction epsilon 
      \param[in] absLimit absolut limit of residual  
      \param[in] maxIter maximal iteration steps 
      \param[in] verbose verbosity 
    */
    CGInverseOp( const OperatorType &op,
                 double  redEps,
                 double absLimit,
                 int maxIter,
                 int verbose )
    : operator_( op ),
      solver_( absLimit, maxIter, (verbose > 0) )
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
  };

} // end namespace Dune

#endif
