#ifndef DUNE_FEM_LDLSOLVER_HH
#define DUNE_FEM_LDLSOLVER_HH

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/matrix/colcompspmatrix.hh>
#include <dune/fem/solver/inverseoperatorinterface.hh>

#if HAVE_SUITESPARSE_LDL

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


template< class DiscreteFunction,
          class Matrix = SparseRowMatrix< typename DiscreteFunction::DiscreteFunctionSpaceType::RangeFieldType > >
class LDLInverseOperator;

template< class DiscreteFunction, class Matrix >
struct LDLInverseOperatorTraits
{
  typedef DiscreteFunction    DiscreteFunctionType;
  typedef AdaptiveDiscreteFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > SolverDiscreteFunctionType;

  typedef Dune::Fem::Operator< DiscreteFunction, DiscreteFunction > OperatorType;
  typedef OperatorType  PreconditionerType;

  typedef Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction, Matrix > AssembledOperatorType;

  typedef LDLInverseOperator< DiscreteFunction, Matrix >  InverseOperatorType;
  typedef SolverParameter SolverParameterType;
};



/** \class LDLInverseOperator
 *  \ingroup DirectSolver
 *  \brief The %LDL direct sparse solver
 *   Details on %LDL can be found on
 *   http://www.cise.ufl.edu/research/sparse/ldl/
 *  \note This will only work if dune-fem has been configured to use LDL
 */
template< class DF, class Matrix >
class LDLInverseOperator : public InverseOperatorInterface< LDLInverseOperatorTraits< DF, Matrix > >
{
  typedef LDLInverseOperatorTraits< DF, Matrix > Traits;
  typedef InverseOperatorInterface< Traits > BaseType;

  friend class InverseOperatorInterface< Traits >;
public:
  /** \brief this solver does not offer preconditioning option */
  static const bool preconditioningAvailable = false;

  typedef LDLInverseOperator< DF, Matrix > ThisType;

  typedef typename BaseType :: SolverDiscreteFunctionType
    SolverDiscreteFunctionType;

  typedef typename BaseType :: OperatorType           OperatorType;
  typedef typename BaseType :: AssembledOperatorType  AssembledOperatorType;

  // \brief The column-compressed matrix type.
  typedef ColCompMatrix<typename AssembledOperatorType::MatrixType::MatrixBaseType,int> CCSMatrixType;

  typedef typename SolverDiscreteFunctionType::DofType DofType;
  typedef typename SolverDiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  /** \brief Constructor.
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   *  \param[in] verbose verbosity
   */
  LDLInverseOperator(const double& redEps, const double& absLimit, const int& maxIter, const bool& verbose,
        const ParameterReader &parameter = Parameter::container() ) :
    verbose_(verbose && Parameter::verbose(Parameter::solverStatistics )), ccsmat_()
  {}

  /** \brief Constructor.
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   */
  LDLInverseOperator(const double& redEps, const double& absLimit, const int& maxIter,
        const ParameterReader &parameter = Parameter::container() ) :
    verbose_(parameter.getValue<bool>("fem.solver.verbose",false) && Parameter::verbose(Parameter::solverStatistics )),
    ccsmat_()
  {}

  LDLInverseOperator(const SolverParameter &parameter = SolverParameter(Parameter::container()) )
  : BaseType(parameter), verbose_(BaseType::verbose())
  {}

  /** \brief Constructor.
   *  \param[in] op Operator to invert
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   *  \param[in] verbose verbosity
   */
  LDLInverseOperator(const OperatorType& op, const double& redEps, const double& absLimit, const int& maxIter, const bool& verbose,
        const ParameterReader &parameter = Parameter::container() ) :
    verbose_(verbose), ccsmat_(), isloaded_(false)
  {
    bind(op);
  }

  /** \brief Constructor.
   *  \param[in] op Operator to invert
   *  \param[in] redEps relative tolerance for residual (not used here)
   *  \param[in] absLimit absolut solving tolerance for residual (not used here)
   *  \param[in] maxIter maximal number of iterations performed (not used here)
   */
  LDLInverseOperator(const OperatorType& op, const double& redEps, const double& absLimit, const int& maxIter,
        const ParameterReader &parameter = Parameter::container() ) :
    verbose_(parameter.getValue<bool>("fem.solver.verbose",false)), ccsmat_(), isloaded_(false)
  {
    bind(op);
  }

  LDLInverseOperator(const OperatorType& op, const ParameterReader &parameter = Parameter::container() ) :
    verbose_(parameter.getValue<bool>("fem.solver.verbose",false)), ccsmat_(), isloaded_(false)
  {
    bind(op);
  }

  // \brief Destructor.
  ~LDLInverseOperator()
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

  // \brief Decompose matrix.
  template<typename... A>
  void prepare(A... ) const
  {
    if( ! assembledOperator_ )
      DUNE_THROW(NotImplemented,"LDLInverseOperator only works for assembled systems!");

    if(!isloaded_)
    {
      ccsmat_ = assembledOperator_->exportMatrix();
      decompose();
      isloaded_ = true;
    }
  }

  // \brief Free allocated memory.
  virtual void finalize()
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
  template<typename... DFs>
  void apply(const TupleDiscreteFunction<DFs...>& arg,TupleDiscreteFunction<DFs...>& dest) const
  {
    // copy DOF's arg into a consecutive vector
    std::vector<DofType> vecArg(arg.size());
    auto vecArgIt(vecArg.begin());
    Hybrid::forEach(std::make_index_sequence<sizeof...(DFs)>{},
      [&](auto i){vecArgIt=std::copy(std::get<i>(arg).dbegin(),std::get<i>(arg).dend(),vecArgIt);});
    std::vector<DofType> vecDest(dest.size());
    // apply operator
    apply(vecArg.data(),vecDest.data());
    // copy back solution into dest
    auto vecDestIt(vecDest.begin());
    Hybrid::forEach(std::make_index_sequence<sizeof...(DFs)>{},[&](auto i){for(auto& dof:dofs(std::get<i>(dest))) dof=(*(vecDestIt++));});
  }

protected:
  void printTexInfo(std::ostream& out) const
  {
    out<<"Solver: LDL direct solver"<<std::endl;
  }

  // \brief Print some statistics about the LDL decomposition.
  void printDecompositionInfo() const
  {
    amd_info(info_);
  }

  /** \brief Get factorization diagonal matrix D.
   *  \warning It is up to the user to preserve consistency.
   */
  DofType* getD()
  {
    return D_;
  }

  /** \brief Get factorization Lp.
   *  \warning It is up to the user to preserve consistency.
   */
  int* getLp()
  {
    return Lp_;
  }

  /** \brief Get factorization Li.
   *  \warning It is up to the user to preserve consistency.
   */
  int* getLi()
  {
    return Li_;
  }

  /** \brief Get factorization Lx.
   *  \warning It is up to the user to preserve consistency.
   */
  DofType* getLx()
  {
    return Lx_;
  }

  /** \brief Get CCS matrix of the operator to solve.
   *  \warning It is up to the user to preserve consistency.
   */
  CCSMatrixType& getCCSMatrix()
  {
    return ccsmat_;
  }

protected:
  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  void apply(const DofType* arg, DofType* dest) const
  {
    prepare();

    // apply part of the call
    const std::size_t dimMat(ccsmat_.N());
    ldl_perm(dimMat, Y_, const_cast<DofType*>(arg), P_);
    ldl_lsolve(dimMat, Y_, Lp_, Li_, Lx_);
    ldl_dsolve(dimMat, Y_, D_);
    ldl_ltsolve(dimMat, Y_, Lp_, Li_, Lx_);
    ldl_permt(dimMat, dest, Y_, P_);

    const_cast<ThisType*>(this)->finalize();
  }

  /** \brief Solve the system.
   *  \param[in] arg right hand side
   *  \param[out] dest solution
   *  \warning You have to decompose the matrix before calling the apply (using the method prepare)
   *   and you have free the decompistion when is not needed anymore (using the method finalize).
   */
  int apply(const SolverDiscreteFunctionType& arg,
             SolverDiscreteFunctionType& dest) const
  {
    apply(arg.leakPointer(),dest.leakPointer());
    return 0;
  }


protected:
  using BaseType :: assembledOperator_;
  const bool verbose_;

  mutable CCSMatrixType ccsmat_;
  mutable bool isloaded_ = false;
  mutable int* Lp_;
  mutable int* Parent_;
  mutable int* Lnz_;
  mutable int* Flag_;
  mutable int* Pattern_;
  mutable int* P_;
  mutable int* Pinv_;
  mutable DofType* D_;
  mutable DofType* Y_;
  mutable DofType* Lx_;
  mutable int* Li_;
  mutable double info_[AMD_INFO];

  // \brief Computes the LDL decomposition.
  void decompose() const
  {
    // allocate vectors
    const std::size_t dimMat(ccsmat_.N());
    D_ = new DofType [dimMat];
    Y_ = new DofType [dimMat];
    Lp_ = new int [dimMat + 1];
    Parent_ = new int [dimMat];
    Lnz_ = new int [dimMat];
    Flag_ = new int [dimMat];
    Pattern_ = new int [dimMat];
    P_ = new int [dimMat];
    Pinv_ = new int [dimMat];

    if(amd_order (dimMat, ccsmat_.getColStart(), ccsmat_.getRowIndex(), P_, (DofType *) NULL, info_) < AMD_OK)
      DUNE_THROW(InvalidStateException,"LDL Error: AMD failed!");
    if(verbose_)
      printDecompositionInfo();
    // compute the symbolic factorisation
    ldl_symbolic(dimMat, ccsmat_.getColStart(), ccsmat_.getRowIndex(), Lp_, Parent_, Lnz_, Flag_, P_, Pinv_);
    // initialise those entries of additionalVectors_ whose dimension is known only now
    Lx_ = new DofType [Lp_[dimMat]];
    Li_ = new int [Lp_[dimMat]];
    // compute the numeric factorisation
    const std::size_t k(ldl_numeric(dimMat, ccsmat_.getColStart(), ccsmat_.getRowIndex(), ccsmat_.getValues(),
                                    Lp_, Parent_, Lnz_, Li_, Lx_, D_, Y_, Pattern_, Flag_, P_, Pinv_));
    // free temporary vectors
    delete [] Flag_;
    delete [] Pattern_;
    delete [] Parent_;
    delete [] Lnz_;

    if(k!=dimMat)
    {
      std::cerr<<"LDL Error: D("<<k<<","<<k<<") is zero!"<<std::endl;
      DUNE_THROW(InvalidStateException,"LDL Error: factorisation failed!");
    }
  }
};

// deprecated old type
template<class DF, class Op, bool symmetric=false>
using LDLOp = LDLInverseOperator< DF, typename Op::MatrixType >;

} // end namespace Fem
} // end namespace Dune

#endif // #if HAVE_SUITESPARSE_LDL

#endif // #ifndef DUNE_FEM_LDLSOLVER_HH
