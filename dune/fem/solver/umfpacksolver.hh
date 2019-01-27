#ifndef DUNE_FEM_UMFPACKSOLVER_HH
#define DUNE_FEM_UMFPACKSOLVER_HH

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/utility.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/colcompspmatrix.hh>

#if HAVE_SUITESPARSE_UMFPACK
#include <dune/fem/misc/umfpack.hh>

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

template<class DiscreteFunction>
class UMFPACKInverseOperator;

template<class DiscreteFunction>
struct UMFPACKInverseOperatorTraits
{
  typedef DiscreteFunction    DiscreteFunctionType;
  typedef AdaptiveDiscreteFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > SolverDiscreteFunctionType;

  typedef Dune::Fem::Operator< DiscreteFunction, DiscreteFunction > OperatorType;
  typedef OperatorType  PreconditionerType;

  typedef Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction > AssembledOperatorType;

  typedef UMFPACKInverseOperator< DiscreteFunction>  InverseOperatorType;
  typedef SolverParameter SolverParameterType;
};

/** \class UMFPACKInverseOperator
 *  \ingroup DirectSolver
 *  \brief The %UMFPack direct sparse solver
 *  %UMFPack will always go double precision and supports complex numbers.
 *  Details on UMFPack can be found on http://www.cise.ufl.edu/research/sparse/umfpack/
 *  \note This will only work if dune-fem has been configured to use UMFPACK
 */
template<class DiscreteFunction>
class UMFPACKInverseOperator : public InverseOperatorInterface< UMFPACKInverseOperatorTraits< DiscreteFunction > >
{
  public:
  typedef UMFPACKInverseOperatorTraits< DiscreteFunction > Traits;
  typedef InverseOperatorInterface< Traits > BaseType;

  typedef typename BaseType :: SolverDiscreteFunctionType
    SolverDiscreteFunctionType;

  typedef typename BaseType :: OperatorType           OperatorType;
  typedef typename BaseType :: AssembledOperatorType  AssembledOperatorType;

  typedef DiscreteFunction DiscreteFunctionType;

  // \brief The column-compressed matrix type.
  typedef ColCompMatrix<typename AssembledOperatorType::MatrixType::MatrixBaseType> CCSMatrixType;
  typedef typename DiscreteFunctionType::DofType DofType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  /** \brief Constructor.
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   *  \param[in] verbose verbosity
   */
  UMFPACKInverseOperator(const double& redEps, const double& absLimit, const int& maxIter, const bool& verbose,
            const ParameterReader &parameter = Parameter::container() ) :
    UMFPACKInverseOperator(verbose)
  {}

  /** \brief Constructor.
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   */
  UMFPACKInverseOperator(const double& redEps, const double& absLimit,
            const int& maxIter, const ParameterReader& parameter = Parameter::container() ) :
    UMFPACKInverseOperator(parameter.getValue<bool>("fem.solver.verbose",false))
  {}
  UMFPACKInverseOperator(const double& redEps, const double& absLimit,
            const ParameterReader& parameter = Parameter::container() ) :
    UMFPACKInverseOperator(parameter.getValue<bool>("fem.solver.verbose",false))
  {}

  UMFPACKInverseOperator(const ParameterReader& parameter ) :
    UMFPACKInverseOperator(parameter.getValue<bool>("fem.solver.verbose",false))
  {}

  UMFPACKInverseOperator(const SolverParameter &parameter = SolverParameter(Parameter::container()) ) :
    verbose_( parameter.verbose() && Dune::Fem::Parameter::verbose() )
  {}

  /** \brief Constructor.
   *  \param[in] op Operator to invert
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   *  \param[in] verbose verbosity
   */
  UMFPACKInverseOperator(const OperatorType& op, const double& redEps, const double& absLimit, const int& maxIter, const bool& verbose,
            const ParameterReader &parameter = Parameter::container() ) :
    UMFPACKInverseOperator(verbose)
  {
    bind(op);
  }

  /** \brief Constructor.
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   */
  UMFPACKInverseOperator(const OperatorType& op, const double& redEps, const double& absLimit,
            const int& maxIter, const ParameterReader& parameter = Parameter::container() ) :
    UMFPACKInverseOperator(parameter.getValue<bool>("fem.solver.verbose",false))
  {
    bind(op);
  }

  UMFPACKInverseOperator(const OperatorType& op, const ParameterReader& parameter = Parameter::container() ) :
    UMFPACKInverseOperator(parameter.getValue<bool>("fem.solver.verbose",false))
  {
    bind(op);
  }

  // \brief Destructor.
  ~UMFPACKInverseOperator()
  {
    unbind();
  }

  void bind( const OperatorType& op )
  {
    // clear old storage
    finalize();
    BaseType::bind( op );
  }

  void unbind ()
  {
    finalize();
    BaseType::unbind();
  }

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
    if(assembledOperator_ && !isloaded_)
    {
      ccsmat_ = std::make_unique<CCSMatrixType>(assembledOperator_->systemMatrix().matrix());
      decompose();
      isloaded_ = true;
    }
  }

  // \brief Free allocated memory.
  void finalize() const
  {
    if(isloaded_)
    {
      getCCSMatrix().free();
      ccsmat_ = nullptr;
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
  int apply(const DofType* arg, DofType* dest) const
  {
    double UMF_Apply_Info[UMFPACK_INFO];
    Caller::solve(UMFPACK_A, getCCSMatrix().getColStart(), getCCSMatrix().getRowIndex(), getCCSMatrix().getValues(),
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
    return 1;
  }

  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  int apply(const AdaptiveDiscreteFunction<DiscreteFunctionSpaceType>& arg,
             AdaptiveDiscreteFunction<DiscreteFunctionSpaceType>& dest) const
  {
    return apply(arg.leakPointer(),dest.leakPointer());
  }

  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  int apply(const ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpaceType>& arg,
             ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpaceType>& dest) const
  {
    // copy DOF's arg into a consecutive vector
    std::vector<DofType> vecArg(arg.size());
    std::copy(arg.dbegin(),arg.dend(),vecArg.begin());
    std::vector<DofType> vecDest(dest.size());
    // apply operator
    apply(vecArg.data(),vecDest.data());
    // copy back solution into dest
    std::copy(vecDest.begin(),vecDest.end(),dest.dbegin());
    return 1;
  }

  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  template<typename... DFs>
  int apply(const TupleDiscreteFunction<DFs...>& arg,TupleDiscreteFunction<DFs...>& dest) const
  {
    // copy DOF's arg into a consecutive vector
    std::vector<DofType> vecArg(arg.size());
    auto vecArgIt(vecArg.begin());
    Hybrid::forEach(Std::make_index_sequence<sizeof...(DFs)>{},
      [&](auto i){vecArgIt=std::copy(std::get<i>(arg).dbegin(),std::get<i>(arg).dend(),vecArgIt);});
    std::vector<DofType> vecDest(dest.size());
    // apply operator
    apply(vecArg.data(),vecDest.data());
    // copy back solution into dest
    auto vecDestIt(vecDest.begin());
    Hybrid::forEach(Std::make_index_sequence<sizeof...(DFs)>{},[&](auto i){for(auto& dof:dofs(std::get<i>(dest))) dof=(*(vecDestIt++));});
    return 1;
  }

  void printTexInfo(std::ostream& out) const
  {
    out<<"Solver: UMFPACK direct solver"<<std::endl;
  }

  double averageCommTime() const
  {
    return 0.0;
  }

  int iterations() const
  {
    return 0;
  }
  void setMaxIterations ( int ) {}

  /** \brief Get CCS matrix of the operator to solve.
   *  \warning It is up to the user to preserve consistency.
   */
  CCSMatrixType& getCCSMatrix() const
  {
    return *ccsmat_;
  }

  // \brief Print some statistics about the UMFPACK decomposition.
  void printDecompositionInfo() const
  {
    Caller::report_info(UMF_Control,UMF_Decomposition_Info);
  }

  using BaseType::assembledOperator_;
  const bool verbose_;
  mutable std::unique_ptr<CCSMatrixType> ccsmat_;
  mutable bool isloaded_ = false;
  mutable void *UMF_Symbolic;
  mutable void *UMF_Numeric;
  mutable double UMF_Control[UMFPACK_CONTROL];
  mutable double UMF_Decomposition_Info[UMFPACK_INFO];

  UMFPACKInverseOperator(const UMFPACKInverseOperator &other)
    : BaseType(other), verbose_(other.verbose_),
      ccsmat_(nullptr),
      UMF_Symbolic(other.UMF_Symbolic),
      UMF_Numeric(other.UMF_Numeric)
  {
    for (int i=0;i<UMFPACK_CONTROL;++i) UMF_Control[i] = other.UMF_Control[i];
    for (int i=0;i<UMFPACK_INFO;++i) UMF_Decomposition_Info[i] = other.UMF_Decomposition_Info[i];
  }

  private:
  explicit UMFPACKInverseOperator ( const bool &verbose ) :
    verbose_(verbose), ccsmat_(nullptr)
  {
    Caller::defaults(UMF_Control);
    UMF_Control[UMFPACK_PRL] = 4;
  }

  typedef typename Dune::UMFPackMethodChooser<DofType> Caller;

  // \brief Computes the UMFPACK decomposition.
  void decompose() const
  {
    const std::size_t dimMat(getCCSMatrix().N());
    Caller::symbolic(static_cast<int>(dimMat), static_cast<int>(dimMat), getCCSMatrix().getColStart(), getCCSMatrix().getRowIndex(),
                     reinterpret_cast<double*>(getCCSMatrix().getValues()), &UMF_Symbolic, UMF_Control, UMF_Decomposition_Info);
    Caller::numeric(getCCSMatrix().getColStart(), getCCSMatrix().getRowIndex(), reinterpret_cast<double*>(getCCSMatrix().getValues()),
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

// deprecated old type
template<class DF, class Op, bool symmetric=false>
using UMFPACKOp = UMFPACKInverseOperator< DF >;

}
}

#endif // #if HAVE_SUITESPARSE_UMFPACK

#endif // #ifndef DUNE_FEM_UMFPACKSOLVER_HH
