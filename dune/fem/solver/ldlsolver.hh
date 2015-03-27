#ifndef DUNE_FEM_LDLSOLVER_HH
#define DUNE_FEM_LDLSOLVER_HH

#include <limits>
#include <iostream>
#include <type_traits>

#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction/blockvectorfunction.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/matrix/colcompspmatrix.hh>

#if HAVE_LDL
#ifdef __cplusplus
extern "C"
{
#include "ldl.h"
#include "amd.h"
}
#endif

namespace Dune
{
namespace Fem
{

/** @addtogroup DirectSolver
 *
 *  In this section implementations of direct solvers
 *  for solving linear systems of the from
 *  \f$A x = b\f$, where \f$A\f$ is a Mapping or
 *  Operator and \f$x\f$ and \f$b\f$ are discrete functions
 *  (see DiscreteFunctionInterface) can be found.
 */

/** \class LDLOp
 *  \ingroup DirectSolver
 *  \brief The %LDL direct sparse solver
 *   Details on %LDL can be found on
 *   http://www.cise.ufl.edu/research/sparse/ldl/
 *  \note This will only work if dune-fem has been configured to use LDL
 */
template<class DF, class Op>
class LDLOp:public Operator<DF, DF>
{
  public:
  typedef DF DiscreteFunctionType;
  typedef Op OperatorType;

  private:
  // \brief The column-compressed matrix type.
  typedef ColCompMatrix<typename OperatorType::MatrixType> CCSMatrixType;
  typedef typename DiscreteFunctionType::DofType DofType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  public:
  /** \brief Constructor.
   *  \param[in] op Operator to invert
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   *  \param[in] verbose verbosity
   */
  LDLOp(const OperatorType& op, const double& redEps, const double& absLimit, const int& maxIter, const bool& verbose) :
    op_(op), verbose_(verbose), ccsmat_(), isloaded_(false)
  {}

  /** \brief Constructor.
   *  \param[in] op Operator to invert
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   */
  LDLOp(const OperatorType& op, const double& redEps=0.0, const double& absLimit=0.0, const int& maxIter=std::numeric_limits<int>::max()) :
    op_(op), verbose_(Parameter::getValue<bool>("fem.solver.verbose",false)), ccsmat_(), isloaded_(false)
  {}

  // \brief Destructor.
  ~LDLOp()
  {}

  /** \brief Solve the system
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   */
  void operator()(const DiscreteFunctionType& arg, DiscreteFunctionType& dest) const
  {
    prepare();
    apply(arg,dest);
    finalize();
  }

  // \brief Decompose matrix.
  template<typename... A>
  void prepare(A... ) const
  {
    if(!isloaded_)
    {
      ccsmat_ = op_.systemMatrix().matrix();
      decompose();
      isloaded_ = true;
    }
  }

  // \brief Free allocated memory.
  inline void finalize() const
  {
    if(isloaded_)
    {
      ccsmat_.free();
      delete [] D_;
      delete [] Y_;
      delete [] Lp_;
      delete [] Lx_;
      delete [] Li_;
      delete [] P_;
      delete [] Pinv_;
      isloaded_ = false;
    }
  }

  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  void apply(const DofType*& arg, DofType*& dest) const
  {
    const std::size_t dimMat(ccsmat_.N());
    ldl_perm(dimMat, Y_, arg, P_);
    ldl_lsolve(dimMat, Y_, Lp_, Li_, Lx_);
    ldl_dsolve(dimMat, Y_, D_);
    ldl_ltsolve(dimMat, Y_, Lp_, Li_, Lx_);
    ldl_permt(dimMat, dest, Y_, P_);
  }

  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  void apply(const AdaptiveDiscreteFunction<DiscreteFunctionSpaceType>& arg,
             AdaptiveDiscreteFunction<DiscreteFunctionSpaceType>& dest) const
  {
    const DofType* argPtr(arg.leakPointer());
    DofType* destPtr(dest.leakPointer());
    apply(argPtr,destPtr);
  }

  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  void apply(const ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpaceType>& arg,
             ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpaceType>& dest) const
  {
    const DofType* argPtr(arg.allocDofPointer());
    DofType* destPtr(dest.allocDofPointer());
    apply(argPtr,destPtr);
    arg.freeDofPointerNoCopy(argPtr);
    dest.freeDofPointer(destPtr);
  }

  inline void printTexInfo(std::ostream& out) const
  {
    out<<"Solver: LDL direct solver";
    out<<"\\\\ \n";
  }

  inline double averageCommTime() const
  {
    return 0.0;
  }

  inline int iterations() const
  {
    return 0;
  }

  /** \brief Get factorization diagonal matrix D.
   *  \warning It is up to the user to preserve consistency.
   */
  inline double*& getD()
  {
    return D_;
  }

  /** \brief Get factorization Lp.
   *  \warning It is up to the user to preserve consistency.
   */
  inline int*& getLp()
  {
    return Lp_;
  }

  /** \brief Get factorization Li.
   *  \warning It is up to the user to preserve consistency.
   */
  inline int*& getLi()
  {
    return Li_;
  }

  /** \brief Get factorization Lx.
   *  \warning It is up to the user to preserve consistency.
   */
  inline double*& getLx()
  {
    return Lx_;
  }

  private:
  const OperatorType& op_;
  const bool verbose_;
  mutable CCSMatrixType ccsmat_;
  mutable bool isloaded_;
  mutable int* Lp_;
  mutable int* Parent_;
  mutable int* Lnz_;
  mutable int* Flag_;
  mutable int* Pattern_;
  mutable int* P_;
  mutable int* Pinv_;
  mutable double* D_;
  mutable double* Y_;
  mutable double* Lx_;
  mutable int* Li_;

  // /brief Computes the LDL decomposition.
  void decompose() const
  {
    // allocate vectors
    const std::size_t dimMat(ccsmat_.N());
    D_ = new double [dimMat];
    Y_ = new double [dimMat];
    Lp_ = new int [dimMat + 1];
    Parent_ = new int [dimMat];
    Lnz_ = new int [dimMat];
    Flag_ = new int [dimMat];
    Pattern_ = new int [dimMat];
    P_ = new int [dimMat];
    Pinv_ = new int [dimMat];

    double Info [AMD_INFO];
    if(amd_order (dimMat, ccsmat_.getColStart(), ccsmat_.getRowIndex(), P_, (double *) NULL, Info) < AMD_OK)
      std::cout<<"WARNING: call to AMD failed."<<std::endl;
    if(verbose_)
      amd_info (Info);
    // compute the symbolic factorisation
    ldl_symbolic(dimMat, ccsmat_.getColStart(), ccsmat_.getRowIndex(), Lp_, Parent_, Lnz_, Flag_, P_, Pinv_);
    // initialise those entries of additionalVectors_ whose dimension is known only now
    Lx_ = new double [Lp_[dimMat]];
    Li_ = new int [Lp_[dimMat]];
    // compute the numeric factorisation
    const int rank(ldl_numeric(dimMat, ccsmat_.getColStart(), ccsmat_.getRowIndex(), ccsmat_.getValues(),
                               Lp_, Parent_, Lnz_, Li_, Lx_, D_, Y_, Pattern_, Flag_, P_, Pinv_));
    // free temporary vectors
    delete [] Flag_;
    delete [] Pattern_;
    delete [] Parent_;
    delete [] Lnz_;

    if(rank!=dimMat)
      std::cout<<"WARNING: matrix is singular."<<std::endl;
  }
};

}
}

#endif // #if HAVE_SPQR

#endif // #ifndef DUNE_FEM_SPQRSOLVER_HH
