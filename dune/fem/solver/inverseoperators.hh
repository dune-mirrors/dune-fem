#ifndef DUNE_FEM_INVERSEOPERATORS_HH
#define DUNE_FEM_INVERSEOPERATORS_HH

#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>

namespace Dune
{

  // CGInverseOperator
  // -----------------
 
  /** \class   CGInverseOperator
   *  \ingroup OEMSolver
   *  \brief   Inverse operator base on CG method
   */
  template< class DiscreteFunction, 
            class Op = Fem::Operator< DiscreteFunction, DiscreteFunction > >
  class CGInverseOperator
  : public Fem::CGInverseOperator< DiscreteFunction >
  {
    typedef Fem::CGInverseOperator< DiscreteFunction > BaseType;
    
  public:
    typedef typename BaseType::DomainFunctionType DomainFunctionType;
    typedef typename BaseType::RangeFunctionType RangeFunctionType;

    typedef DomainFunctionType DestinationType;
    
    //! type of operator 
    typedef Op OperatorType;
    
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
                        double redEps, double absLimit,
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
                        double redEps, double absLimit,
                        unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
    : BaseType( op, redEps, absLimit, maxIter ),
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
    CGInverseOperator ( const OperatorType &op,
                        const PreconditioningType &precond,
                        double redEps, double absLimit,
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
      if( preconditioning && LinearOperator :: assembled ) 
      {
        // create diagonal preconditioner 
        precondObj_ = new Fem::DiagonalPreconditioner< DomainFunctionType, LinearOperator >( linearOp );
        preconditioner_ = precondObj_;
      }
    }

    using BaseType::preconditioner_;
    PreconditioningType *precondObj_;
  };

#ifdef DUNE_FEM_COMPATIBILITY
// to be removed after release of version 1.3 
#define CGInverseOp  CGInverseOperator
#endif

} // namespace Dune

#endif // #ifndef DUNE_FEM_INVERSEOPERATORS_HH
