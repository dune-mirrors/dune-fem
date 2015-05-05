#ifndef DUNE_FEM_UMFPACKSOLVER_HH
#define DUNE_FEM_UMFPACKSOLVER_HH

#include <limits>
#include <iostream>
#include <type_traits>
#include <vector>

#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction/blockvectorfunction.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/matrix/colcompspmatrix.hh>


#if HAVE_DUNE_ISTL
#include <dune/istl/umfpack.hh>

#ifdef ENABLE_UMFPACK

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

/** \class SPQROp
 *  \ingroup DirectSolver
 *  \brief The %UMFPack direct sparse solver
 *  %UMFPack will always go double precision and supports complex numbers.
 *  Details on UMFPack can be found on http://www.cise.ufl.edu/research/sparse/umfpack/
 *  \note This will only work if dune-fem has been configured to use UMFPACK
 */
template<class DF, class Op, bool symmetric=false>
class UMFPACKOp:public Operator<DF, DF>
{
  public:
  typedef DF DiscreteFunctionType;
  typedef Op OperatorType;

  // \brief The column-compressed matrix type.
  typedef ColCompMatrix<typename OperatorType::MatrixType::MatrixBaseType> CCSMatrixType;
  typedef typename DiscreteFunctionType::DofType DofType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  /** \brief Constructor.
   *  \param[in] op Operator to invert
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   *  \param[in] verbose verbosity
   */
  UMFPACKOp(const OperatorType& op, const double& redEps, const double& absLimit, const int& maxIter, const bool& verbose) :
    op_(op), verbose_(verbose), ccsmat_(), isloaded_(false)
  {
    Caller::defaults(UMF_Control);
    UMF_Control[UMFPACK_PRL] = 4;
  }

  /** \brief Constructor.
   *  \param[in] op Operator to invert
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   */
  UMFPACKOp(const OperatorType& op, const double& redEps=0.0, const double& absLimit=0.0,
            const int& maxIter=std::numeric_limits<int>::max()) :
    op_(op), verbose_(Parameter::getValue<bool>("fem.solver.verbose",false)), ccsmat_(), isloaded_(false)
  {
    Caller::defaults(UMF_Control);
    UMF_Control[UMFPACK_PRL] = 4;
  }

  // \brief Destructor.
  ~UMFPACKOp()
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
      Caller::free_symbolic(&UMF_Symbolic);
      Caller::free_numeric(&UMF_Numeric);
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
    double UMF_Apply_Info[UMFPACK_INFO];
    Caller::solve(UMFPACK_A, ccsmat_.getColStart(), ccsmat_.getRowIndex(), ccsmat_.getValues(),
                  dest, const_cast<DofType*>(arg), UMF_Numeric, UMF_Control, UMF_Apply_Info);
    if(verbose_)
    {
      Caller::report_status(UMF_Control, UMF_Apply_Info[UMFPACK_STATUS]);
      std::cout <<"[UMFPack Solve]" << std::endl;
      std::cout << "Wallclock Time: " << UMF_Apply_Info[UMFPACK_SOLVE_WALLTIME]
                << " (CPU Time: " << UMF_Apply_Info[UMFPACK_SOLVE_TIME] << ")" << std::endl;
      std::cout << "Flops Taken: " << UMF_Apply_Info[UMFPACK_SOLVE_FLOPS] << std::endl;
      std::cout << "Iterative Refinement steps taken: " << UMF_Apply_Info[UMFPACK_IR_TAKEN] << std::endl;
      std::cout << "Error Estimate: " << UMF_Apply_Info[UMFPACK_OMEGA1] << " resp. " << UMF_Apply_Info[UMFPACK_OMEGA2] << std::endl;
    }
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
    // convert ISTLBlockVectorDiscreteFunction to AdaptiveDiscreteFunction in order to have a DOF*
    std::vector<DofType> vecArg(arg.size());
    AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> adaptiveArg(arg.name(),arg.space(),vecArg.data());
    adaptiveArg.assign(arg);
    std::vector<DofType> vecDest(dest.size());
    AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> adaptiveDest(dest.name(),dest.space(),vecDest.data());
    adaptiveDest.assign(dest);

    apply(adaptiveArg,adaptiveDest);

    // copy solution into dest
    dest.assign(adaptiveDest);
  }

  inline void printTexInfo(std::ostream& out) const
  {
    out<<"Solver: UMFPACK direct solver";
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

  /** \brief Get CCS matrix of the operator to solve.
   *  \warning It is up to the user to preserve consistency.
   */
  inline CCSMatrixType& getCCSMatrix()
  {
    return ccsmat_;
  }

  // /brief Print some statistics about the UMFPACK decomposition.
  inline void printDecompositionInfo() const
  {
    Caller::report_info(UMF_Control,UMF_Decomposition_Info);
  }

  private:
  const OperatorType& op_;
  const bool verbose_;
  mutable CCSMatrixType ccsmat_;
  mutable bool isloaded_;
  mutable void *UMF_Symbolic;
  mutable void *UMF_Numeric;
  mutable double UMF_Control[UMFPACK_CONTROL];
  mutable double UMF_Decomposition_Info[UMFPACK_INFO];

  typedef typename Dune::UMFPackMethodChooser<DofType> Caller;

  // /brief Computes the UMFPACK decomposition.
  void decompose() const
  {
    const std::size_t dimMat(ccsmat_.N());
    Caller::symbolic(static_cast<int>(dimMat), static_cast<int>(dimMat), ccsmat_.getColStart(), ccsmat_.getRowIndex(),
                     reinterpret_cast<double*>(ccsmat_.getValues()), &UMF_Symbolic, UMF_Control, UMF_Decomposition_Info);
    Caller::numeric(ccsmat_.getColStart(), ccsmat_.getRowIndex(), reinterpret_cast<double*>(ccsmat_.getValues()),
                    UMF_Symbolic, &UMF_Numeric, UMF_Control, UMF_Decomposition_Info);
    if(verbose_)
    {
      Caller::report_status(UMF_Control,UMF_Decomposition_Info[UMFPACK_STATUS]);
      std::cout << "[UMFPack Decomposition]" << std::endl;
      std::cout << "Wallclock Time taken: " << UMF_Decomposition_Info[UMFPACK_NUMERIC_WALLTIME]
                << " (CPU Time: " << UMF_Decomposition_Info[UMFPACK_NUMERIC_TIME] << ")" << std::endl;
      std::cout << "Flops taken: " << UMF_Decomposition_Info[UMFPACK_FLOPS] << std::endl;
      std::cout << "Peak Memory Usage: " << UMF_Decomposition_Info[UMFPACK_PEAK_MEMORY]*UMF_Decomposition_Info[UMFPACK_SIZE_OF_UNIT]
                << " bytes" << std::endl;
      std::cout << "Condition number estimate: " << 1./UMF_Decomposition_Info[UMFPACK_RCOND] << std::endl;
      std::cout << "Numbers of non-zeroes in decomposition: L: " << UMF_Decomposition_Info[UMFPACK_LNZ]
                << " U: " << UMF_Decomposition_Info[UMFPACK_UNZ] << std::endl;
    }
  }
};

}
}

#endif // #if HAVE_DUNE_ISTL

#endif // #if HAVE_UMFPACK

#endif // #ifndef DUNE_FEM_UMFPACKSOLVER_HH
