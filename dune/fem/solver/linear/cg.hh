#ifndef DUNE_FEM_CG_HH
#define DUNE_FEM_CG_HH

#include <type_traits>
#include <vector>
#include <cmath>

#include <dune/common/ftraits.hh>
#include <dune/fem/solver/parameter.hh>

namespace Dune
{
namespace Fem
{
namespace LinearSolver
{

  template <class Operator, class Precoditioner, class DiscreteFunction>
  inline int cg( Operator &op,
                 Precoditioner* preconditioner,
                 std::vector< DiscreteFunction >& tempMem,
                 DiscreteFunction& x,
                 const DiscreteFunction& b,
                 const double epsilon,
                 const int maxIterations,
                 const int toleranceCriteria,
                 std::ostream* os = nullptr )
  {
    typedef typename DiscreteFunction::RangeFieldType RangeFieldType;
    typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

    assert( preconditioner ? tempMem.size() == 5 : tempMem.size() == 3 );

    DiscreteFunction& h = tempMem[ 0 ];

    //h=Ax
    op( x, h );

    //r=Ax-b
    DiscreteFunction& r = tempMem[ 1 ];
    r.assign( h );
    r -= b;

    //p=b-A*x <= r_0 Deufelhard
    DiscreteFunction& p = tempMem[ 2 ];
    p.assign( b );
    p -= h;

    DiscreteFunction& s = ( preconditioner ) ? tempMem[ 3 ] : p;

    //q=B*p <=q Deufelhard
    DiscreteFunction& q = ( preconditioner ) ? tempMem[ 4 ] : p;
    if( preconditioner )
    {
      (*preconditioner)( p, q );
      s.assign( q );
    }

    RangeFieldType prevResidual = 0;    // note that these will be real_type but require scalar product evaluation
    RangeFieldType residual = p.scalarProductDofs( q );//<p,Bp>

    const RealType tolerance = epsilon * epsilon * (
      toleranceCriteria == ToleranceCriteria::relative ? b.normSquaredDofs( ) :
      toleranceCriteria == ToleranceCriteria::residualReduction ? std::real(residual) :
      // absolute tolerance
      1.0
    );

    int iterations = 0;
    for( iterations = 0; (std::real(residual) > tolerance) && (iterations < maxIterations); ++iterations )
    {
      if( iterations > 0 )
      {
        assert( residual/prevResidual == residual/prevResidual );
        const RangeFieldType beta= residual / prevResidual;
        q *= beta;
        if( preconditioner )
        {
          q += s;
        }
        else
        {
          p -= r;
        }
      }

      op( q, h );

      RangeFieldType qdoth = q.scalarProductDofs( h );
      const RangeFieldType alpha = residual / qdoth;//<p,Bp>/<q,Aq>
      assert( !std::isnan( alpha ) );
      x.axpy( alpha, q );

      if( preconditioner )
      {
        p.axpy( -alpha, h );//r_k
        (*preconditioner)( p, s); //B*r_k

        prevResidual = residual;//<rk-1,B*rk-1>
        residual = p.scalarProductDofs( s );//<rk,B*rk>
      }
      else
      {
        r.axpy( alpha, h );

        prevResidual = residual;
        residual = r.normSquaredDofs( );
      }

      if( os )
      {
        (*os) << "Fem::CG it: " << iterations << " : residual " << std::sqrt(residual) << std::endl;
      }
    }

    return (iterations < maxIterations) ? iterations : -iterations;
  }

} // naemspace LinearSolver
} // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_CGINVERSEOPERATOR_HH
