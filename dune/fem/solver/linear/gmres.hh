#ifndef DUNE_FEM_GMRES_HH
#define DUNE_FEM_GMRES_HH

#include <cmath>
#include <cassert>
#include <iostream>

#include <utility>

#include <dune/fem/solver/parameter.hh>

namespace Dune
{
namespace Fem
{
namespace LinearSolver
{

  namespace detail {
    template <class T>
    class Matrix
    {
    public:
      Matrix(int n, int m) : n_(n), m_(m), data_( n*m, T(0) ) {}
      Matrix(int n) : Matrix(n,n) {}

      // element access
      T& operator()(int i, int j)
      {
        assert(i>=0 && i<n_ && j>=0 && j<m_);
        return data_[i*m_ + j];
      }

      // const element access
      T operator()(int i, int j) const
      {
        assert(i>=0 && i<n_ && j>=0 && j<m_);
        return data_[i*m_ + j];
      }

      // conversion operators
      operator T*() { return data_.data(); }
      operator const T *() const { return data_.data(); }

    protected:
      const int n_, m_;
      std::vector< T > data_;
    };
  }

  //! return x * y
  template <class FieldType>
  FieldType scalarProduct( const int dim, const FieldType *x, const FieldType* y )
  {
    FieldType scp = 0;
    for( int i=0; i<dim; ++i )
    {
      scp += x[ i ] * y[ i ];
    }
    return scp;
  }

  // computes y = beta y + alpha op(A) x
  template <class Communication, class FieldType, class DiscreteFunction>
  void gemv(const Communication& comm,
            const int m,           // j+1
            std::vector< DiscreteFunction >& v,
            const DiscreteFunction& vjp,
            FieldType *y           // global_dot
           )
  {
      for(int l=0; l<m; ++l)
      {
        y[ l ] = 0;
      }

      const auto& auxiliaryDofs = vjp.space().auxiliaryDofs();
      const auto& vj = vjp.dofVector();

      auto scp = [&m, &y, &vj, &v] (const size_t dof)
      {
        for(int l=0; l<m; ++l)
        {
          y[ l ] += (vj[ dof ] * v[ l ].dofVector()[ dof ]);
        }
      };
      // see dune/fem/space/common/auiliarydofs.hh
      forEachPrimaryDof( auxiliaryDofs, scp );

      // communicate sum
      comm.sum( y, m );
  }

  //! dblas_rotate with inc = 1
  template<class FieldType>
  void rotate( const int dim,
               FieldType* x, FieldType* y,
               const FieldType c, const FieldType s)
  {
    int i = dim;
    while (i--)
    {
      const FieldType _x = *x;
      const FieldType _y = *y;
      *x = c*_x + s*_y;
      *y = c*_y - s*_x;
      ++x;
      ++y;
    }
  }

  // Saad, Youcef;  Schultz, Martin H.
  // GMRES: A generalized minimal residual algorithm for solving nonsymmetric
  // linear systems. (English)
  // [J] SIAM J. Sci. Stat. Comput. 7, 856-869 (1986). [ISSN 0196-5204]
  template <class Operator, class Preconditioner, class DiscreteFunction>
  inline int gmres( Operator& op, Preconditioner* preconditioner,
                    std::vector< DiscreteFunction >& v,
                    DiscreteFunction& u,
                    const DiscreteFunction& b,
                    const int m, // gmres inner iterations
                    const double tolerance,
                    const int maxIterations,
                    const int toleranceCriteria,
                    std::ostream* os = nullptr )
  {
    typedef typename DiscreteFunction :: RangeFieldType FieldType;

    const auto& comm = u.space().gridPart().comm();

    detail::Matrix< FieldType > H( m+1, m ); // \in \R^{m+1 \times m}
    std::vector< FieldType > g_( 6*m, 0.0 );

    FieldType* g = g_.data();
    FieldType* s = g + (m+1);
    FieldType* c = s + m;
    FieldType* y = c + m;

    DiscreteFunction& v0 = v[ 0 ];

    std::vector< FieldType > global_dot( m+1, FieldType(0) );

    // relative or absolute tolerance
    double _tolerance = tolerance;
    if (toleranceCriteria == ToleranceCriteria::relative)
    {
      global_dot[ 0 ] = b.scalarProductDofs( b );
      _tolerance *= std::sqrt(global_dot[0]);
    }

    int iterations = 0;
    while (true)
    {
      // start
      op(u, v0);

      v0 -= b ;

      // scalarProduct( dim, v, v );
      // contains 1 allreduce (global zum)
      global_dot[ 0 ] = v0.scalarProductDofs( v0 );

      //comm.allreduce(1, local_dot, global_dot, MPI_SUM);
      FieldType res = std::sqrt(global_dot[0]);

      if (toleranceCriteria == ToleranceCriteria::residualReduction && iterations==0)
      {
        _tolerance *= res;
      }

      if (os)
      {
        (*os) << "Fem::GMRES outer iteration : " << res << std::endl;
      }

      if (res < _tolerance) break;

      g[0] = -res;
      for(int i=1; i<=m; i++) g[i] = 0.0;

      // cblas_dscal(dim, 1.0/res, v, 1);
      v0 *= (1.0/res);

      //scale( dim, 1.0/res, v );

      // iterate
      for(int j=0; j<m; j++)
      {
        DiscreteFunction& vj  = v[ j ];
        DiscreteFunction& vjp = v[ j + 1 ];

        // apply the linear operator (perhaps in combination with the
        // preconditioner)
        if (preconditioner)
        {
          DiscreteFunction& z = v[ m+1 ];
          (*preconditioner)(vj, z );
          op( z, vjp);
        }
        else
        {
          op(vj, vjp);
        }

        // contains 1 allreduce (global sum)
        gemv(comm, j+1, v, vjp, global_dot.data());

        for(int i=0; i<=j; i++) H(i,j) = global_dot[i];

        // assuming beta == 1.0
        for(int l=0; l<j+1; ++l)
        {
          vjp.axpy( -global_dot[l], v[l] );
        }

        // 1 allreduce (global sum)
        global_dot[ 0 ] = vjp.scalarProductDofs( vjp );

        H(j+1,j) = std::sqrt(global_dot[0]);
        // scale(dim, 1.0/H(j+1,j), vjp );
        vjp *= 1.0/H(j+1,j);

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

        if ( os )
        {
          (*os) << "Fem::GMRES it: " << iterations << " : " <<  std::abs(g[j+1]) << std::endl;
        }

        ++iterations;
        if (std::abs(g[j+1]) < _tolerance
            || iterations >= maxIterations ) break;
      }

      //
      // form the approximate solution
      //

      int last = iterations%m;
      if (last == 0) last = m;

      // compute y via backsubstitution
      for(int i=last-1; i>=0; --i)
      {
        const FieldType dot = scalarProduct( last-(i+1), &H(i,i)+1, &y[i+1] );
        y[i] = (g[i] - dot)/ H(i,i);
      }

      // update the approx. solution
      if (preconditioner)
      {
        // u += M^{-1} (v[0], ..., v[last-1]) y
        DiscreteFunction& u_tmp = v[ m ]; // we don't need this vector anymore
        DiscreteFunction& z = v[ m+1 ];
        u_tmp.clear();

        // u += (v[0], ..., v[last-1]) y
        for(int i=0; i<last; ++i)
        {
          u_tmp.axpy( y[ i ], v[ i ] );
        }

        (*preconditioner)(u_tmp, z);
        u += z;
      }
      else{
        // u += (v[0], ..., v[last-1]) y
        for(int i=0; i<last; ++i)
        {
          u.axpy( y[ i ], v[ i ] );
        }
      }

      if (std::abs(g[last]) < _tolerance) break;
    }

    // output
    if ( os ) {
      (*os) << "Fem::GMRES: number of iterations: "
         << iterations
         << std::endl;
    }

    return (iterations < maxIterations) ? iterations : -iterations;
  }


} // end namespace Solver

} // end namespace Fem

} // end namespace Dune

#endif
