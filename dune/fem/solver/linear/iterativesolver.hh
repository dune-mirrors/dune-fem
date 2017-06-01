#ifndef DUNE_FEM_ITERATIVE_SOLVER_HH
#define DUNE_FEM_ITERATIVE_SOLVER_HH

#include <cmath>
#include <cassert>
#include <iostream>

#include <utility>

#include <dune/fem/storage/dynamicarray.hh>

namespace Dune
{

namespace Fem
{

namespace LinearSolver
{

  struct ToleranceCriteria
  {
    static const int absolute = 0;
    static const int relative = 1;
    static const int residualReduction = 2;
  };

  /** \brief set tolerance for linear solvers
   *  \param solver      linear solver object
   *  \param redEps      relative reduction tolerance
   *  \param absLimit    absolute error tolerance
   *  \param paramName   parameter key for Fem::Parameter
   */
  template <class Solver>
  inline
  void setTolerance(const Dune::Fem::ParameterReader &parameter, Solver &solver,
                    double redEps, double absLimit, const char *paramName)
  {
    const std::string errorTypeTable[] =
      { "absolute", "relative", "residualreduction" };
    const int errorType = parameter.getEnum( paramName, errorTypeTable, 0 );
    solver.setTolerance( absLimit, errorType );
  }

  template <class T>
  class FunctionIF
  {
  protected:
    FunctionIF() {}
  public:
    virtual void operator()( const T* x, T* y ) = 0;
    virtual int size() const = 0 ;
    virtual ~FunctionIF() {}
  };

  template <class T, class Comm >
  class IterativeSolver
  {
    IterativeSolver(const IterativeSolver &) = delete;
  public:
    typedef T         FieldType;
    typedef FieldType field_type;

    typedef Fem::StaticArray< FieldType > VectorType;
    typedef Fem::DynamicArray< FieldType > MutableVectorType;

    typedef FunctionIF< FieldType > Function;

    IterativeSolver(const Comm &comm) :
      comm_( comm ),
      preconditioner_( nullptr ),
      os_( nullptr ),
      tolerance( 1e-8 ),
      max_num_of_iterations( 100 ),
      num_of_iterations( 0 ),
      toleranceCriteria( ToleranceCriteria::absolute )
    {
    }

    virtual int number_of_iterations() { return num_of_iterations; }

    virtual void setTolerance( const double tol, const int crit )
    {
      tolerance = tol;
      toleranceCriteria = crit;
    }

    virtual void set_max_number_of_iterations( const int max_it )
    {
      max_num_of_iterations = max_it;
    }

    virtual void set_output( std::ostream& os )
    {
      os_ = &os ;
    }

    // from Function, solve Au = b, Au = op(u)
    virtual bool solve(Function &op, FieldType *u, const FieldType *b) = 0;

    virtual void set_preconditioner(Function &preconditioner)
    {
      preconditioner_ = &preconditioner;
    }

    // set pointer to preconditioner to zero
    virtual void unset_preconditioner()
    {
      preconditioner_ = nullptr;
      z_.clear();
    }

  protected:
    //! return x * y
    FieldType scalarProduct( const int dim, const FieldType *x, const FieldType* y ) const
    {
      FieldType scp = 0;
      for( int i=0; i<dim; ++i )
      {
        scp += x[ i ] * y[ i ];
      }
      return comm_.sum( scp );
    }

    //! scale x with factor alpha
    void scale( const int dim, const FieldType alpha, FieldType *x ) const
    {
      for( int i=0; i<dim; ++i )
      {
        x[ i ] *= alpha;
      }
    }

    //! dblas_rotate with inc = 1
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

    // op(A) \in \R^{m \times n}
    // x \in \R^n
    // y \in \R^m
    //
    // computes y = beta y + alpha op(A) x
    void gemv(const int m, const int n,
              const FieldType alpha, const FieldType *A, const int lda,
              const FieldType *x,
              // const int incx,
              const FieldType beta, FieldType *y
              // const int incy
             ) const
    {
      assert(lda >= n); // consistent lda

      for(int i=0; i<m; ++i)
      {
        double sum = 0.0;
        for(int j=0; j<n; ++j)
        {
          sum += A[lda*i + j] * x[j];
        }
        y[i] = beta * y[i] + alpha * sum;
      }

      // communicate sum
      comm_.sum( y, m );
    }


    // A \in \R^{m \times k}, lda>=k
    // B \in \R^{k \times n}, ldb>=n
    // C \in \R^{m \times n}, ldc>=n
    //
    // computes C = alpha A*B + beta C
    void gemm(
        //const enum CBLAS_ORDER order,
        // const enum CBLAS_TRANSPOSE transA,
        // const enum CBLAS_TRANSPOSE transB,
         const int m, const int n, const int k,
         const FieldType alpha, const FieldType *A, const int lda,
         const FieldType *B, const int ldb,
         const FieldType beta, FieldType *C, const int ldc) const
    {
      // new
      // consistent lda, ldb, ldc
      assert(lda >= k && ldb >= n && ldc >= n);

      for(int i=0; i<m; i++)
      {
        for(int j=0; j<n; j++)
        {
          double sum = 0.0;
          for(int l=0; l<k; l++)
          {
            sum += A[i*lda + l] * B[l*ldb + j];
          }
          C[i*ldc + j] = alpha*sum + beta*C[i*ldc + j];
        }
      }
    }

    virtual void resize(int new_size)
    {
      if (preconditioner_)
      {
        z_.resize( new_size, 0.0 );
      }
    }

  protected:
    const Comm& comm_;
    Function* preconditioner_;
    std::ostream* os_;

    double tolerance;

    int max_num_of_iterations ;
    int num_of_iterations ;
    int toleranceCriteria ;

    MutableVectorType z_;
  };

} // end namespace Solver

} // end namespace Fem

} // end namespace Dune

#endif
