// GMRES inline implementation

#ifndef DUNE_FEM_GMRES_SIMPLE_HH
#define DUNE_FEM_GMRES_SIMPLE_HH

#include <cmath>
#include <cassert>
#include <iostream>

#include <utility>

#include <dune/fem/solver/linear/iterativesolver.hh>
#include <dune/fem/solver/pardg.hh>

namespace Dune
{
namespace Fem
{
namespace LinearSolver
{
  // Saad, Youcef;  Schultz, Martin H.
  // GMRES: A generalized minimal residual algorithm for solving nonsymmetric
  // linear systems. (English)
  // [J] SIAM J. Sci. Stat. Comput. 7, 856-869 (1986). [ISSN 0196-5204]
  template <class T, class Comm >
  class GMRES : public IterativeSolver< T, Comm >
  {
  protected:
    typedef IterativeSolver< T, Comm >  BaseType;
    using BaseType :: scalarProduct;
    using BaseType :: os_;
    using BaseType :: preconditioner_;
    using BaseType :: tolerance;
    using BaseType :: max_num_of_iterations;
    using BaseType :: num_of_iterations;
    using BaseType :: toleranceCriteria;
    using BaseType :: scale;
    using BaseType :: rotate;
    using BaseType :: gemv;
    using BaseType :: gemm;
    using BaseType :: z_;

    GMRES(const GMRES &) = delete;

  public:
    typedef typename BaseType :: Function  Function;

    typedef T         FieldType;
    typedef FieldType field_type;

    typedef typename BaseType :: VectorType         VectorType;
    typedef typename BaseType :: MutableVectorType  MutableVectorType;

    GMRES(const Comm &comm, const int _m) :
      BaseType( comm ),
      m(_m),
      H(m+1,m),
      g_( 6*m, 0.0 ),
      v_()
    {
      g_.setMemoryFactor( 1.1 );
      v_.setMemoryFactor( 1.1 );
    }

    GMRES(GMRES &&other) :
      BaseType( other ),
      m(other.m),
      H(std::move(other.H)),
      g_( std::move( other.g_ )),
      v_( std::move(other.v_))
    {
    }

    // from Function, solve Au = b, Au = op(u)
    bool solve(Function &op, FieldType* u, const FieldType* b);

  protected:
    virtual void resize(int new_size)
    {
      BaseType :: resize( new_size );
      v_.resize( (m+1)*new_size, 0.0 );
    }

  protected:
    const int m;
    PARDG::Matrix H; // \in \R^{m+1 \times m}
    MutableVectorType g_, v_;
  };

  template <class T, class Comm>
  inline bool GMRES< T, Comm>::solve(Function &op, FieldType* u, const FieldType *b)
  {
    const int dim = op.size();
    resize(dim);

    FieldType* g = g_.data();
    FieldType* s = g + (m+1);
    FieldType* c = s + m;
    FieldType* y = c + m;
    FieldType* global_dot = y + m;

    FieldType* v = v_.data();

    // relative or absolute tolerance
    double _tolerance = tolerance;
    if (toleranceCriteria == ToleranceCriteria::relative){
      //local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
      //comm.allreduce(1, local_dot, global_dot, MPI_SUM);

      global_dot[ 0 ] = scalarProduct( dim, b, b );
      _tolerance *= std::sqrt(global_dot[0]);
    }

    FieldType* z = ( preconditioner_ ) ? z_.data() : nullptr;

    int iterations = 0;
    while (true)
    {
      // start
      op(u, v);
      // cblas_daxpy(dim, -1.0, b, 1, v, 1);
      for( int i=0; i<dim; ++i )
        v[ i ] -= b[ i ];

      global_dot[ 0 ] = scalarProduct( dim, v, v );

      //comm.allreduce(1, local_dot, global_dot, MPI_SUM);
      FieldType res = std::sqrt(global_dot[0]);

      if (toleranceCriteria == ToleranceCriteria::residualReduction && iterations==0)
      {
        _tolerance *= res;
      }

      if (os_)
      {
        (*os_) << "GMRES outer iteration : " << res << std::endl;
      }

      if (res < _tolerance) break;

      g[0] = -res;
      for(int i=1; i<=m; i++) g[i] = 0.0;

      // cblas_dscal(dim, 1.0/res, v, 1);
      scale( dim, 1.0/res, v );

      // iterate
      for(int j=0; j<m; j++){
        FieldType *vj = v + j*dim;
        FieldType *vjp = vj + dim;

        // apply the linear operator (perhaps in combination with the
        // preconditioner)
        if (preconditioner_)
        {
          (*preconditioner_)(vj, z);
          op(z, vjp);
        }
        else op(vj, vjp);

        //cblas_dgemv(CblasRowMajor, CblasNoTrans,
        //            j+1, dim, 1.0, v, dim, vjp, 1, 0.0, global_dot, 1);
                    //j+1, dim, 1.0, v, dim, vjp, 1, 0.0, local_dot, 1);
        gemv(j+1, dim, 1.0, v, dim, vjp, 0.0, global_dot);
        //comm.allreduce(j+1, local_dot, global_dot, MPI_SUM);

        for(int i=0; i<=j; i++) H(i,j) = global_dot[i];

        //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        //            1, dim, j+1,  -1.0, global_dot, m,  v, dim,  1.0, vjp, dim);
        gemm(1, dim, j+1,  -1.0, global_dot, m,  v, dim,  1.0, vjp, dim);

        //local_dot[0] = cblas_ddot(dim, vjp, 1, vjp, 1);
        //comm.allreduce(1, local_dot, global_dot, MPI_SUM);
        global_dot[ 0 ] = scalarProduct( dim, vjp, vjp );

        H(j+1,j) = std::sqrt(global_dot[0]);
        // cblas_dscal(dim, 1.0/H(j+1,j), vjp, 1);
        scale(dim, 1.0/H(j+1,j), vjp );

        // perform Givens rotation
        for(int i=0; i<j; i++)
        {
          rotate(1, &H(i+1,j), &H(i,j), c[i], s[i]);
        }
        const FieldType h_j_j = H(j,j);
        const FieldType h_jp_j = H(j+1,j);
        const FieldType norm = std::sqrt(h_j_j*h_j_j + h_jp_j*h_jp_j);
        c[j] = h_j_j / norm;
        s[j] = -h_jp_j / norm;
        rotate(1, &H(j+1,j), &H(j,j), c[j], s[j]);
        rotate(1, &g[j+1], &g[j], c[j], s[j]);

        if ( os_ )
        {
          (*os_) << "GMRES it: " << iterations << " : " <<  std::abs(g[j+1]) << std::endl;
        }

        iterations++;
        if (std::abs(g[j+1]) < _tolerance
            || iterations >= max_num_of_iterations) break;
      }

      //
      // form the approximate solution
      //

      int last = iterations%m;
      if (last == 0) last = m;
      // compute y via backsubstitution
      for(int i=last-1; i>=0; i--)
      {
        const FieldType dot = scalarProduct( last-(i+1), &H(i,i)+1, &y[i+1] );
        y[i] = (g[i] - dot)/ H(i,i);
      }

      // update the approx. solution
      if (preconditioner_)
      {
        // u += M^{-1} (v[0], ..., v[last-1]) y
        FieldType* u_tmp = v + m*dim; // we don't need this vector anymore
        //dset(dim, 0.0, u_tmp, 1);
        for( int i=0; i<dim; ++i ) u_tmp[ i ] = 0.0;

        for(int i=0; i<last; i++)
        {
          FieldType *vi = v + i*dim;
          const FieldType alpha = y[i];
          // cblas_daxpy(dim, y[i], vi, 1, u_tmp, 1);
          for( int k=0; k<dim; ++k )
            u_tmp[ k ] += alpha * vi[ k ];
        }
        (*preconditioner_)(u_tmp, z);
        // cblas_daxpy(dim, 1.0, z, 1, u, 1);
        for( int i=0; i<dim; ++i )
          u[ i ] += z[ i ];
      }
      else{
        // u += (v[0], ..., v[last-1]) y
        for(int i=0; i<last; i++)
        {
          FieldType *vi = v + i*dim;
          const FieldType alpha = y[i];
          // cblas_daxpy(dim, y[i], vi, 1, u, 1);
          for( int k=0; k<dim; ++k )
            u[ k ] += alpha * vi[ k ];
        }
      }

      if (std::abs(g[last]) < _tolerance) break;
    }

    // output
    if ( os_ ) {
      (*os_) << "GMRES: number of iterations: "
         << iterations
         << std::endl;
    }

    // update the global number of iterations from IterativeSolver
    num_of_iterations += iterations;

    return (iterations < max_num_of_iterations)? true: false;
  }


} // end namespace Solver

} // end namespace Fem

} // end namespace Dune

#endif
