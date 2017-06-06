#ifndef DUNE_FEM_BICGSTAB_HH
#define DUNE_FEM_BICGSTAB_HH

#include <cmath>
#include <cassert>
#include <iostream>

#include <utility>

#include <dune/fem/solver/linear/iterativesolver.hh>

namespace Dune
{
namespace Fem
{
namespace LinearSolver
{

  template <class Communication, class DiscreteFunction, class FieldType>
  void scalarProductVecs( const Communication& comm,
                          const DiscreteFunction& r_df,
                          const DiscreteFunction& s_df,
                          const DiscreteFunction& r_star_df,
                          FieldType* global_dot )
  {
    const auto& r      = r_df.dofVector();
    const auto& s      = s_df.dofVector();
    const auto& r_star = r_star_df.dofVector();

    global_dot[ 4 ] = global_dot[ 3 ] = global_dot[ 2 ] = global_dot[ 1 ] = global_dot[ 0 ] = 0.0;

    const size_t size = r.size();
    for( size_t i=0; i<size; ++i )
    {
      global_dot[0] += r[ i ] * s[ i ];
      global_dot[1] += r[ i ] * r[ i ];
      global_dot[2] += s[ i ] * s[ i ];
      global_dot[3] += s[ i ] * r_star[ i ];
      global_dot[4] += r[ i ] * r_star[ i ];
    }

    // make sum global
    comm.sum( global_dot, 5 );
  }

  template <class Communication, class DofVector>
  typename DofVector::FieldType
  scalarProduct( const Communication& comm,
                 const DofVector& x,
                 const DofVector& y )
  {
    typename DofVector::FieldType scp = 0;

    const size_t size = x.size();
    for( size_t i=0; i<size; ++i )
    {
      scp += x[ i ] * y[ i ];
    }

    // make sum global
    return comm.sum( scp );
  }

  /* General implementation of a BiCG-stab algorithm based on Dune::Fem::DiscreteFunction
   * \param op linear operator A
   * \param x  solution to be sought
   * \param b  right hand side
   * \param tempMem temporary memory
   * \param tolerance for the solver
   * \param maxIter maximal iterations to be performed
   * \param toleranceCriteria tolerance citeria (see iterativesolvers.hh)
   * \param os output if enabled
   */
  template <class Operator, class Precoditioner, class DiscreteFunction>
  inline int bicgstab( Operator &op,
                       Precoditioner* preconditioner,
                       std::vector< DiscreteFunction >& tempMem,
                       DiscreteFunction& x,
                       const DiscreteFunction& b,
                       const double tolerance,
                       const int maxIterations,
                       const int toleranceCriteria,
                       std::ostream* os = nullptr )
  {
    const auto& comm = b.space().gridPart().comm();

    DiscreteFunction& r = tempMem[ 0 ];
    DiscreteFunction& r_star = tempMem[ 1 ];
    DiscreteFunction& p = tempMem[ 2 ];
    DiscreteFunction& s = tempMem[ 3 ];
    DiscreteFunction& tmp = tempMem[ 4 ];

    DiscreteFunction* z = nullptr;

    typedef typename DiscreteFunction :: RangeFieldType FieldType;

    if( tempMem.size() == 6 )
    {
      z = &tempMem[ 5 ];
    }

    // relative or absolute tolerance
    FieldType global_dot[5];
    double _tolerance = tolerance;

    if (toleranceCriteria == ToleranceCriteria::relative)
    {
      global_dot[ 0 ] = scalarProduct( comm, b.dofVector(), b.dofVector() );
      _tolerance *= std::sqrt(global_dot[0]);
    }

    // init
    if (preconditioner)         // right preconditioning
    {
      (*preconditioner)(x, *z); // z = M^{-1} x
      op( *z, r);               // r = A z
    }
    else
    {
      op(x, r);              // r = A x
    }

    // residual: r = b - r
    r *= -1.0;
    r += b ;

    // p = r
    p.assign( r );
    // r_star = r
    r_star.assign( r );

    global_dot[0] = scalarProduct( comm, r.dofVector(), r_star.dofVector() );

    FieldType nu = global_dot[0];
    if (toleranceCriteria == ToleranceCriteria::residualReduction)
    {
      _tolerance *= std::sqrt(nu);
    }

    // iterate
    int iterations = 0;
    while (true)
    {
      // 2x linear operator 1x dot
      if( preconditioner ) // right preconditioning
      {
        (*preconditioner)(p, *z); // z = M^{-1} p
        op( *z, tmp);              // tmp = A z
      }
      else op(p, tmp);           // tmp = A p

      global_dot[0] = scalarProduct( comm, tmp.dofVector(), r_star.dofVector() );

      const FieldType alpha = nu / global_dot[0];
      s.assign( r );
      s.axpy( -alpha, tmp );

      if( preconditioner )
      {
        // right preconditioning
        (*preconditioner)(s, *z); // z = M^{-1} s
        op(*z, r);                 // r = A z
      }
      else
      {
        op(s, r);                 // r = A s
      }

      // 5x dot: r * s | r * r | s * s | s * r_star | r * r_star
      scalarProductVecs( comm, r, s, r_star, global_dot );

      // scalars
      const FieldType omega = global_dot[0] / global_dot[1];
      const FieldType res = std::sqrt(global_dot[2]
            -omega*(2.0*global_dot[0] - omega*global_dot[1]) );

      const FieldType beta = (global_dot[3] - omega*global_dot[4])
        *alpha / (omega*nu);

      nu = (global_dot[3] - omega*global_dot[4]);

      // update
      // x += alpha * p + omega s
      x.axpy( alpha, p );
      x.axpy( omega, s );

      // r = s - omega r
      r *= -omega;
      r += s ;

      // p = r + beta * ( p - omega * tmp )
      p *= beta;
      p.axpy( -omega*beta, tmp );
      p += r;

      ++iterations;
      if (res < _tolerance || iterations >= maxIterations ) break;
      if ( os )
      {
        (*os) << "Fem::BiCGstab it: " << iterations << " : " << res << std::endl;
      }
    }

    // output
    if ( os )
    {
      (*os) << "Fem::BiCGstab:  number of iterations: "
         << iterations
         << std::endl;
    }

    // setup approx solution for right preconditioner use
    if( preconditioner )
    {
      // right preconditioner
      (*preconditioner)(x, *z);

      // x = z
      x.assign( *z );
    }

    if( iterations >= maxIterations )
      return -iterations;

    return iterations;
  }

} // end namespace Solver

} // end namespace Fem

} // end namespace Dune

#endif
