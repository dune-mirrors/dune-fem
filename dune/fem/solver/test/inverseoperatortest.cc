#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// C++ includes
#include <iostream>
#include <string>
#include <type_traits>

// dune-common includes
#include <dune/common/ftraits.hh>

// dune-grid includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

// dune-fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/solver/viennacl.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#if HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#include <dune/fem/solver/istl.hh>
#endif // HAVE_DUNE_ISTL

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscinverseoperators.hh>
#endif // HAVE_PETSC

#if HAVE_EIGEN
#include <dune/fem/storage/eigenvector.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/operator/linear/eigenoperator.hh>
#include <dune/fem/solver/eigen.hh>
#endif // HAVE_EIGEN

#if HAVE_SUITESPARSE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#endif // HAVE_SUITESPARSE_UMFPACK

#if HAVE_SUITESPARSE_LDL
#include <dune/fem/solver/ldlsolver.hh>
#endif // HAVE_SUITESPARSE_LDL

#if HAVE_SUITESPARSE_SPQR
#include <dune/fem/solver/spqrsolver.hh>
#endif // HAVE_SUITESPARSE_SPQR

#include <dune/fem/solver/amgxsolver.hh>

// local includes
#include "../../test/massoperator.hh"

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

static constexpr int dim = 2;
static constexpr int polOrder = 2;

using GridType      = Dune::YaspGrid< dim >;
using GridPartType  = Dune::Fem::LeafGridPart< GridType >;

using FieldType = typename GridType::ctype;
using RealType  = typename Dune::FieldTraits< FieldType >::real_type;

using SpaceType         = Dune::Fem::FunctionSpace< FieldType, FieldType, dim, 1 >;
using DiscreteSpaceType = Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, polOrder >;

struct Function
{
  using FunctionSpaceType = SpaceType;

  using DomainType        = typename FunctionSpaceType::DomainType;
  using RangeType         = typename FunctionSpaceType::RangeType;
  using JacobianRangeType = typename FunctionSpaceType::JacobianRangeType;

  void evaluate ( const DomainType &x, RangeType &value ) const
  {
    value = 1.;
    for( int k = 0; k < FunctionSpaceType::dimDomain; ++k )
      value *= sin( M_PI * x[k] );
  }
};


template< class InverseOperator, class LinearOperator = typename InverseOperator::OperatorType >
struct Algorithm
{
  using InverseOperatorType = InverseOperator;
  using LinearOperatorType  = LinearOperator;

  using DiscreteFunctionType   = typename InverseOperatorType::DomainFunctionType;
  using MassOperatorType       = MassOperator< DiscreteFunctionType, LinearOperatorType >;
  using AffineMassOperatorType = AffineMassOperator< DiscreteFunctionType, LinearOperatorType >;

  static_assert( std::is_base_of< typename InverseOperator::OperatorType, LinearOperator >::value, "type mismatch in Algorithm." );
  static_assert( std::is_same< typename InverseOperator::RangeFunctionType, DiscreteFunctionType >::value, "type mismatch in Algorithm." );

  template <class SolverParam = Dune::Fem::SolverParameter>
  static bool apply( GridType& grid, const std::string& designation, const bool verboseSolver = false,
                     SolverParam *param = nullptr)
  {
    const double eps = 5e-6;

    GridPartType gridPart( grid );
    DiscreteSpaceType space( gridPart );

    DiscreteFunctionType u( "u", space );
    DiscreteFunctionType rhs( "rhs", space );
    u.clear();

    Function f;
    auto gridFunction = Dune::Fem::gridFunctionAdapter( f, gridPart, space.order()+1 );
    MassOperatorType massOperator( space );
    massOperator.assembleRHS(gridFunction, rhs);
    AffineMassOperatorType affineMassOperator( space, gridFunction );

    unsigned long maxIter = 2*space.size();
    maxIter = space.gridPart().comm().sum( maxIter );

    auto f_ = gridFunctionAdapter( "exact", f, gridPart, polOrder+2 );
    Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );

    bool pass = true;

    if (param)
    {
      InverseOperatorType inverseOperator ( *param );
      inverseOperator.bind( massOperator );
      inverseOperator( rhs, u );
      auto dist = l2norm.distance( f_, u );
      pass &= dist < eps;
      if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
        std::cout << designation << " InvOp(myParam)\n" << dist << "\n" << std::endl;

      rhs.clear();
      Dune::Fem::NewtonInverseOperator<LinearOperatorType,InverseOperatorType> newtonInvOp( *param );
      newtonInvOp.bind( affineMassOperator );
      newtonInvOp( rhs, u );
      dist = l2norm.distance( f_, u );
      pass &= dist < eps;
      if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
        std::cout << designation << " NewtonInvOp(myParam)\n" << dist << "\n" << std::endl;
    }
    else
    {
      InverseOperatorType inverseOperator;
      inverseOperator.bind( massOperator );
      inverseOperator( rhs, u );
      auto dist = l2norm.distance( f_, u );
      pass &= dist < eps;
      if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
        std::cout << designation << " inverseOp()\n" << dist << "\n" << std::endl;

      auto param = Dune::Fem::parameterDict("fem.solver", "method",Dune::Fem::SolverParameter::cg);
      InverseOperatorType inverseOperatorA( param );
      inverseOperatorA.bind( massOperator );
      inverseOperatorA( rhs, u );
      dist = l2norm.distance( f_, u );
      pass &= dist < eps;
      if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
        std::cout << designation << " inverseOp(NewSParam)\n" << dist << "\n" << std::endl;

      rhs.clear();
      Dune::Fem::NewtonInverseOperator<LinearOperatorType,InverseOperatorType> newtonInvOp;
      newtonInvOp.bind( affineMassOperator );
      newtonInvOp( rhs, u );
      dist = l2norm.distance( f_, u );
      pass &= dist < eps;
      if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
        std::cout << designation << " NewtonInvOp()\n" << dist << "\n" << std::endl;

      Dune::Fem::NewtonInverseOperator<LinearOperatorType,InverseOperatorType> newtonInvOpA( param );
      newtonInvOpA.bind( affineMassOperator );
      newtonInvOpA( rhs, u );
      dist = l2norm.distance( f_, u );
      pass &= dist < eps;
      if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
        std::cout << designation << " NewtonInvOp(NewSParam)\n" << dist << "\n" << std::endl;

      Dune::Fem::NewtonInverseOperator<LinearOperatorType,InverseOperatorType>
        newtonInvOpB( Dune::Fem::parameterDict( "fem.solver.newton.", "linear.method","cg" ) );
      newtonInvOpB.bind( affineMassOperator );
      newtonInvOpB( rhs, u );
      dist = l2norm.distance( f_, u );
      pass &= dist < eps;
      if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
        std::cout << designation << " NewtonInvOp(paramDict)\n" << dist << "\n" << std::endl;
    }
    return pass;
  }
};

std::shared_ptr< GridType > createGrid( const double length, const int cells )
{
  std::stringstream file;
  file << "DGF" << std::endl;
  file << "Interval" << std::endl;
  for( int i=0; i<dim; ++i )
    file << "0 ";
  file << std::endl;
  for( int i=0; i<dim; ++i )
    file << length << " ";
  file << std::endl;
  for( int i=0; i<dim; ++i )
    file << cells << " ";
  file << std::endl;

  file << "#" << std::endl;
  file << "GridParameter" << std::endl;
  file << "overlap 0" << std::endl;
  file << "#" << std::endl;

  std::shared_ptr< GridType > ptr;
  ptr.reset( Dune::GridPtr< GridType > (file).release() );
  return ptr;
}

int main(int argc, char** argv)
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

  //GridType grid({1., 1.}, {16, 16});
  auto gridPtr = createGrid( 1., 16 );
  GridType& grid = *gridPtr;

  const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();
  grid.globalRefine( step );

  bool verboseSolver = Dune::Fem::SolverParameter().verbose();
  bool pass = true;

  // CGInverseOperator + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::CGInverseOperator< DiscreteFunction >;

    std::string designation(" === CGInverseOperator + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // KrylovInverseOperator + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;

    using InverseOperator   = Dune::Fem::KrylovInverseOperator< DiscreteFunction >;
    std::string designation1(" === KrylovInverseOperator + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation1, verboseSolver );

    designation1 = std::string(" === KrylovInverseOperator + SparseRowLinearOperator + SolverParameter === ");
    Dune::Fem::SolverParameter param( Dune::Fem::parameterDict(
            "fem.solver.",
            "method","cg", "newton.linear.method","gmres"));
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation1, verboseSolver, &param);

    using CgInverseOperator = Dune::Fem::CgInverseOperator< DiscreteFunction >;
    std::string designation2(" === CgInverseOperator + SparseRowLinearOperator === ");
    pass &= Algorithm< CgInverseOperator, LinearOperator >::apply( grid, designation2, verboseSolver );

    using GmresInverseOperator = Dune::Fem::GmresInverseOperator< DiscreteFunction >;
    std::string designation3(" === GmresInverseOperator + SparseRowLinearOperator === ");
    pass &= Algorithm< GmresInverseOperator, LinearOperator >::apply( grid, designation3, verboseSolver );

    using BicgstabInverseOperator = Dune::Fem::BicgstabInverseOperator< DiscreteFunction >;
    std::string designation4(" === BicgstabInverseOperator + SparseRowLinearOperator === ");
    pass &= Algorithm< BicgstabInverseOperator, LinearOperator >::apply( grid, designation4, verboseSolver );
  }

#if HAVE_SUITESPARSE_UMFPACK
  // CGInverseOperator + SparseRowLinearOperator
  if( Dune::Fem::MPIManager::size() == 1 )
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::UMFPACKInverseOperator< DiscreteFunction >;

    std::string designation(" === UMFPACK + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif // HAVE_SUITESPARSE_LDL

#if HAVE_SUITESPARSE_LDL
  // CGInverseOperator + SparseRowLinearOperator
  if( Dune::Fem::MPIManager::size() == 1 )
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::LDLOp< DiscreteFunction, LinearOperator >;

    std::string designation(" === LDLOp + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif // HAVE_SUITESPARSE_LDL

#if 0 // HAVE_SUITESPARSE_SPQR // this fails on the dune CI runner with Debian11 - deactivate for now
  // CGInverseOperator + SparseRowLinearOperator
  if( Dune::Fem::MPIManager::size() == 1 )
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::SPQRInverseOperator< DiscreteFunction >;

    std::string designation(" === SPQROp + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif // HAVE_SUITESPARSE_SPQR

#if HAVE_DUNE_ISTL
  // ISTLInverseOperator + ISTLLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::ISTLLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction >;

    std::string designation(" === ISTLInverseOperator + ISTLLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );

    designation = std::string(" === ISTLInverseOperator + ISTLLinearOperator + ISTLSolverParameter === ");
    Dune::Fem::ISTLSolverParameter param(std::string("istlparam."));
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver,
                    &param );

  }

  // ISTLInverseOperator< ISTLCGSolver > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::SolverParameter::cg >;

    std::string designation(" === ISTLInverseOperator< ISTLCGSolver > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // ISTLInverseOperator< ISTLBiCGSTABSolver > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::SolverParameter::bicgstab >;

    std::string designation(" === ISTLInverseOperator< ISTLBiCGSTABSolver > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // ISTLInverseOperator< ISTLMINRESSolver > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::SolverParameter::minres >;

    std::string designation(" === ISTLInverseOperator< ISTLMINRESSolver > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // ISTLInverseOperator< ISTLRestartedGMRes > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::SolverParameter::gmres >;

    std::string designation(" === ISTLInverseOperator< ISTLRestartedGMRes > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // ISTL::InverseOperator< LinearOperator > + ISTLLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::ISTLLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTL::InverseOperator< LinearOperator >;

    std::string designation(" === ISTL::InverseOperator + ISTLLinearOperator === ");
    InverseOperator::SolverParameterType param( Dune::Fem::parameterDict(
        "istl.",
            "verbosity", verboseSolver ? "full" : "off"
      ) );
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver, &param );
  }

#if HAVE_SUPERLU
  // ISTLInverseOperator< ISTLSuperLU > + ISTLLinearOperator
  if( Dune::Fem::MPIManager::size() == 1 )
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::ISTLLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::SolverParameter::superlu >;

    std::string designation(" === ISTLInverseOperator< ISTLSuperLU > + ISTLLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif // HAVE_SUPERLU

#endif // HAVE_DUNE_ISTL

#if HAVE_PETSC
  // PetscInverseOperator + PetscLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::PetscDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::PetscInverseOperator< DiscreteFunction >;

    std::string designation(" === PetscInverseOperator + PetscLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );

    // REMARK: note that in this version of passing in the parameters both
    // a direct call to the inverse operator and the to the newton inverse
    // operator will use the default prefix of the SolverParameter (i.e.  fem.solver.)
    // In this example the prefix is set to "petsctest" (again for linear and non linear)
    // An empty prefix leads to issues since 'method' without prefix leads to a Warning
    // The same issue happens if one passes  "petsctest." as first argument to parameterDict
    designation = std::string(" === PetscInverseOperator + PetscLinearOperator + PetscParameter === ");
    Dune::Fem::PetscSolverParameter param( "petsctest.", Dune::Fem::parameterDict(
            "petsctest.",
              "preconditioning.method","hypre",
              "petsc.hypre.method", "pilu-t",
              "verbose",false
            ));
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver, &param);
  }

  /*
  // PetscInverseOperator + PetscLinearOperator + AdaptiveDiscreteFunction
  // this case is appearing for adaptive simulations with PETSc as solver backend
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::PetscInverseOperator< DiscreteFunction >;

    std::string designation(" === PetscInverseOperator + PetscLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );

    // REMARK: note that in this version of passing in the parameters both
    // a direct call to the inverse operator and the to the newton inverse
    // operator will use the default prefix of the SolverParameter (i.e.  fem.solver.)
    // In this example the prefix is set to "petsctest" (again for linear and non linear)
    // An empty prefix leads to issues since 'method' without prefix leads to a Warning
    // The same issue happens if one passes  "petsctest." as first argument to parameterDict
    designation = std::string(" === PetscInverseOperator + PetscLinearOperator + PetscParameter === ");
    Dune::Fem::PetscSolverParameter param( "petsctest.", Dune::Fem::parameterDict(
            "petsctest.",
              "preconditioning.method","hypre",
              "petsc.hypre.method", "pilu-t",
              "verbose",false
            ));
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver, &param);
  }
  */
#endif // HAVE_PETSC

  /*
#if HAVE_EIGEN
  // EigenCGInverseOperator + EigenLinearOperator
  {
    using DofVector         = Dune::Fem::EigenVector< FieldType >;
    using DiscreteFunction  = Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DiscreteSpaceType, DofVector > >;
    using LinearOperator    = Dune::Fem::EigenLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::EigenCGInverseOperator< DiscreteFunction >;

    std::string designation(" === EigenCGInverseOperator + EigenLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // EigenBiCGStabInverseOperator + EigenLinearOperator
  {
    using DofVector         = Dune::Fem::EigenVector< FieldType >;
    using DiscreteFunction  = Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DiscreteSpaceType, DofVector > >;
    using LinearOperator    = Dune::Fem::EigenLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::EigenBiCGStabInverseOperator< DiscreteFunction >;

    std::string designation(" === EigenBiCGStabInverseOperator + EigenLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif //HAVE_EIGEN
*/

  /*
#if HAVE_VIENNACL
  // ViennaCLInverseOperator + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ViennaCLInverseOperator< DiscreteFunction >;

    std::string designation(" === ViennaCLInverseOperator + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif
  */

#if HAVE_AMGXSOLVER && HAVE_PETSC
  // AMGX solver wrapper + PetscLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::PetscDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::AMGXInverseOperator< DiscreteFunction >;

    std::string designation(" === AMGXInverseOperator + PetscLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif

  //Dune::Fem::Parameter::write( "param", ".txt", true );

  return pass ? 0 : 1;
}
