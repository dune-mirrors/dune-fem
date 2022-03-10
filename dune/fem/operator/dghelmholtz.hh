#ifndef DUNE_FEM_OPERATOR_DGHELMHOLTZ_HH
#define DUNE_FEM_OPERATOR_DGHELMHOLTZ_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fem/solver/parameter.hh>


namespace Dune
{

  namespace Fem
  {

    // DGHelmholtzJacobianOperator
    // ---------------------------

    template< class JacobianOp >
    class DGHelmholtzJacobianOperator
    : public JacobianOp
    {
      typedef JacobianOp BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

      DGHelmholtzJacobianOperator ( const std::string &name, const DomainFunctionSpaceType &dSpace, const RangeFunctionSpaceType &rSpace,
                                    const SolverParameter& param = SolverParameter() )
      : BaseType( name, dSpace, rSpace ),
        lambda_( 0 ),
        wTmp_( "DGHelmholtzJacobianOperator temporary", rSpace )
      {}

      void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        w.assign( u );
        if( lambda() != 0.0 )
        {
          BaseType::operator()( u, wTmp_ );
          w.axpy( -lambda(), wTmp_ );
        }
      }

      const double &lambda () const { return lambda_; }
      void setLambda ( double lambda ) { lambda_ = lambda; }

    protected:
      double lambda_;
      mutable RangeFunctionType wTmp_;
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
      typedef SpaceOperator SpaceOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::JacobianOperatorType JacobianOperatorType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      explicit DGHelmholtzOperator ( SpaceOperatorType &spaceOp )
      : spaceOp_( spaceOp ),
        lambda_( 0 ),
        wTmp_( "DGHelmholtz::tmp", space() )
      {}

      void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        w.assign( u );
        if( lambda() != 0.0 )
        {
          spaceOperator()( u, wTmp_ );
          w.axpy( -lambda(), wTmp_ );
        }
      }

      void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
      {
        spaceOperator().jacobian( u, jOp );
        jOp.setLambda( lambda() );
      }

      const double &lambda () const { return lambda_; }
      void setLambda ( double lambda ) { lambda_ = lambda; }

      void setTime ( double time ) { spaceOperator().setTime( time ); }

      const DiscreteFunctionSpaceType &space () const { return spaceOperator().space(); }

      void initializeTimeStepSize ( const DomainFunctionType &u ) const
      {
        spaceOperator()( u, wTmp_ );
      }

      double timeStepEstimate () const { return spaceOperator().timeStepEstimate(); }

      const SpaceOperatorType &spaceOperator () const { return spaceOp_; }
      SpaceOperatorType &spaceOperator () { return spaceOp_; }

    protected:
      SpaceOperator &spaceOp_;
      double lambda_;
      mutable RangeFunctionType wTmp_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_DGHELMHOLTZ_HH
