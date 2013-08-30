#ifndef DUNE_FEM_OPERATOR_DGHELMHOLTZ_HH
#define DUNE_FEM_OPERATOR_DGHELMHOLTZ_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>


namespace Dune
{

  namespace Fem
  {

    // DGHelmholtzJacobianOperator
    // ---------------------------

    template< class JacobianOp >
    class DGHelmholtzJacobianOperator
    : public Operator< typename JacobianOp::DomainFunctionType, typename JacobianOp::RangeFunctionType >
    {
      typedef DGHelmholtzJacobianOperator< JacobianOp > ThisType;
      typedef Operator< typename JacobianOp::DomainFunctionType, typename JacobianOp::RangeFunctionType > BaseType;

      template< class > friend class DGHelmholtzOperator;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

      DGHelmholtzJacobianOperator ( const std::string &name, const DomainFunctionSpaceType &dSpace, const RangeFunctionSpaceType &rSpace )
      : jacobianOp_( name, dSpace, rSpace ),
        lambda_( 0 ),
        wTmp_( "DGHelmholtzJacobianOperator temporary", rSpace )
      {}

      void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        w.assign( u );
        if( lambda_ != 0.0 )
        {
          spaceOp_( u, wTmp_ );
          w.axpy( -lambda_, wTmp_ );
        }
      }

      void setLambda ( double lambda ) { lambda_ = lambda; }

    protected:
      JacobianOp jacobianOp_;
      double lambda_;
      RangeFunctionType wTmp_;
    };



    // DGHelmholtzOperator
    // -------------------

    template< class SpaceOperator >
    class DGHelmholtzOperator
    : public DifferentiableOperator< DGHelmholtzJacobianOperator< typename SpaceOperator::JacobianOperatorType > >
    {
      typedef DGHelmholtzOperator< SpaceOperator > ThisType;
      typedef DifferentiableOperator< DGHelmholtzJacobianOperator< typename SpaceOperator::JacobianOperatorType > > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::JacobianOperatorType JacobianOperatorType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      explicit DGHelmholtzOperator ( SpaceOperator &spaceOp )
      : spaceOp_( spaceOp ),
        lambda_( 0 ),
        wTmp_( "DGHelmholtzOperator temporary", space() )
      {}

      void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        w.assign( u );
        if( lambda_ != 0.0 )
        {
          spaceOp_( u, wTmp_ );
          w.axpy( -lambda_, wTmp_ );
        }
      }

      void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
      {
        spaceOp_.jacobian( u, jOp.jacobianOp_ );
        jOp.setLambda( lambda_ );
      }

      void setLambda ( double lambda ) { lambda_ = lambda; }
      void setTime ( double time ) { spaceOp_.setTime( time ); }

      const DiscreteFunctionSpaceType &space () const { return spaceOp_.space(); }

      double timeStepEstimate () const { return spaceOp_.timeStepEstimate(); }

    protected:
      SpaceOperator spaceOp_;
      double lambda_;
      RangeFunctionType wTmp_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_DGHELMHOLTZ_HH
