// BICGSTAB inline implementation

#ifndef DUNE_FEM_BICGSTAB_SIMPLE_HH
#define DUNE_FEM_BICGSTAB_SIMPLE_HH

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

  // Iterative methods for sparse linear systems / Yousef Saad.
  // Boston, Mass. : PWS Publ., 1996. - XVI, 447 S. : graph. Darst.;
  // (engl.)(The PWS series in computer science)
  // ISBN 0-534-94776-X
  //
  // BICGSTAB - Solver page 220
  template <class T, class Comm >
  class BICGSTAB : public IterativeSolver< T, Comm >
  {
  protected:
    typedef IterativeSolver< T, Comm >  BaseType;

    using BaseType :: scalarProduct;
    using BaseType :: os_;
    using BaseType :: preconditioner_;
    using BaseType :: comm_;
    using BaseType :: tolerance;
    using BaseType :: max_num_of_iterations;
    using BaseType :: num_of_iterations;
    using BaseType :: toleranceCriteria;
    using BaseType :: z_;

    BICGSTAB(const BICGSTAB &) = delete;

  public:
    typedef typename BaseType :: Function  Function;

    typedef T         FieldType;
    typedef FieldType field_type;

    typedef typename BaseType :: VectorType         VectorType;
    typedef typename BaseType :: MutableVectorType  MutableVectorType;

    BICGSTAB(const Comm &comm ) :
      BaseType( comm ),
      v_()
    {
      v_.setMemoryFactor( 1.1 );
    }

    BICGSTAB(BICGSTAB &&other) :
      BaseType( other ),
      v_( std::move(other.v_))
    {
    }

    // from Function, solve Au = b, Au = op(u)
    bool solve(Function &op, FieldType *u, const FieldType *b);

  protected:
    virtual void resize(int new_size)
    {
      BaseType :: resize( new_size );
      v_.resize( 5*new_size, 0.0 );
    }

  protected:
    MutableVectorType  v_;
  };

  template <class T, class Comm>
  inline bool BICGSTAB< T, Comm>::solve(Function &op, FieldType* x, const FieldType *b)
  {
    const int dim = op.size();
    resize(dim);

    FieldType* r = v_.data();
    FieldType* r_star = r + dim;
    FieldType* p = r_star + dim;
    FieldType* s = p + dim ;
    FieldType* tmp = s + dim;
    FieldType* z = ( preconditioner_ ) ? z_.data() : nullptr;

    // relative or absolute tolerance
    FieldType global_dot[5];
    double _tolerance = tolerance;
    if (toleranceCriteria == ToleranceCriteria::relative){
      //local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
      //comm.allreduce(1, local_dot, global_dot, MPI_SUM);
      global_dot[ 0 ] = scalarProduct( dim, b, b );
      _tolerance *= sqrt(global_dot[0]);
    }

    // init
    if (preconditioner_){       // right preconditioning
      (*preconditioner_)(x, z); // z = M^{-1} x
      op(z, r);                // r = A z
    }
    else op(x, r);             // r = A x

    for(int k=0; k<dim; k++)
    {
      r[k] = b[k] - r[k];
      p[k] = r[k];
      r_star[k] = r[k];
    }

    //local_dot[0] = cblas_ddot(dim, r, 1, r_star, 1);
    //comm.allreduce(1, local_dot, global_dot, MPI_SUM);
    global_dot[0] = scalarProduct( dim, r, r_star );

    FieldType nu = global_dot[0];
    if (toleranceCriteria == ToleranceCriteria::residualReduction){
      _tolerance *= sqrt(nu);
    }

    // iterate
    int iterations = 0;
    while (true){
      // 2x linear operator 1x dot
      if (preconditioner_){    // right preconditioning
        (*preconditioner_)(p, z); // z = M^{-1} p
        op(z, tmp);              // tmp = A z
      }
      else op(p, tmp);           // tmp = A p

      //local_dot[0] = cblas_ddot(dim, tmp, 1, r_star, 1);
      //comm.allreduce(1, local_dot, global_dot, MPI_SUM);
      global_dot[0] = scalarProduct( dim, tmp, r_star );

      const FieldType alpha = nu / global_dot[0];
      for(int k=0; k<dim; k++) s[k] = r[k] - alpha*tmp[k];

      if (preconditioner_)
      {
        // right preconditioning
        (*preconditioner_)(s, z); // z = M^{-1} s
        op(z, r);                 // r = A z
      }
      else
      {
        op(s, r);                 // r = A s
      }

      // 5x dot
      global_dot[0]=global_dot[1]=global_dot[2]=global_dot[3]=global_dot[4] = 0.0;
      for(int k=0; k<dim; k++){
        global_dot[0] += r[k]*s[k];
        global_dot[1] += r[k]*r[k];
        global_dot[2] += s[k]*s[k];
        global_dot[3] += s[k]*r_star[k];
        global_dot[4] += r[k]*r_star[k];
      }

      //comm.allreduce(5, local_dot, global_dot, MPI_SUM);
      comm_.sum( global_dot, 5 );

      // scalars
      const FieldType omega = global_dot[0] / global_dot[1];
      const FieldType res = sqrt(global_dot[2]
            -omega*(2.0*global_dot[0] - omega*global_dot[1]) );

      const double beta = (global_dot[3] - omega*global_dot[4])
        *alpha / (omega*nu);

      nu = (global_dot[3] - omega*global_dot[4]);

      // update
      for(int k=0; k<dim; k++){
        x[k] += alpha*p[k] + omega*s[k];
        r[k] = s[k] - omega*r[k];
        p[k] = r[k] + beta*( p[k] - omega*tmp[k] );
      }

      iterations++;
      if (res < _tolerance || iterations >= max_num_of_iterations) break;
      if ( os_ )
      {
        (*os_) << "BiCGstab it: " << iterations << " : " << res << std::endl;
      }
    }

    // output
    if ( os_ )
    {
      (*os_) << "BiCGstab:  number of iterations: "
         << iterations
         << std::endl;
    }

    // setup approx solution for right preconditioner use
    if (preconditioner_)
    {
      // right preconditioner
      (*preconditioner_)(x, z);
      // x = z
      for( int i=0; i<dim; ++i )
        x[ i ] = z[ i ];
    }

    // update the global number of iterations from IterativeSolver
    num_of_iterations += iterations;

    return (iterations < max_num_of_iterations)? true: false;
  }


} // end namespace Solver

} // end namespace Fem

} // end namespace Dune

#endif
