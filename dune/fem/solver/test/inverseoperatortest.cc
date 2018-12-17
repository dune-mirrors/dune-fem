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

#if HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
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

#if HAVE_AMGXSOLVER
#include <dune/fem/solver/amgxsolver.hh>
#endif

// local includes
#include "../../test/massoperator.hh"



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

  using DiscreteFunctionType  = typename InverseOperatorType::DomainFunctionType;
  using MassOperatorType      = MassOperator< DiscreteFunctionType, LinearOperatorType >;

  static_assert( std::is_base_of< typename InverseOperator::OperatorType, LinearOperator >::value, "type mismatch in Algorithm." );
  static_assert( std::is_same< typename InverseOperator::RangeFunctionType, DiscreteFunctionType >::value, "type mismatch in Algorithm." );

  static bool apply( GridType& grid, const std::string& designation, bool verboseSolver = false )
  {
    GridPartType gridPart( grid );
    DiscreteSpaceType space( gridPart );

    MassOperatorType massOperator( space );

    DiscreteFunctionType u( "u", space );
    DiscreteFunctionType rhs( "rhs", space );
    u.clear();

    Function f;
    auto gridFunction = Dune::Fem::gridFunctionAdapter( f, gridPart, space.order()+1 );
    massOperator.assembleRHS( gridFunction, rhs );

    unsigned long maxIter = space.size();
    maxIter = space.gridPart().comm().sum( maxIter );

    InverseOperatorType inverseOperator ( 1e-10, 1e-10, maxIter, verboseSolver );
    inverseOperator.bind( massOperator );
    inverseOperator( rhs, u );

    auto f_ = gridFunctionAdapter( "exact", f, gridPart, polOrder+2 );

    Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );

    auto dist = l2norm.distance( f_, u );
    bool pass = dist < 3e-5;

    if( Dune::Fem::Parameter::verbose() || (Dune::Fem::MPIManager::rank() == 0 && !pass) )
      std::cout << designation << "\n" << dist << "\n" << std::endl;

    return pass;
  }
};


int main(int argc, char** argv)
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

  GridType grid({1., 1.}, {4, 4});
  const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();
  grid.globalRefine( 2*step );

  bool verboseSolver = false;
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

#if HAVE_SUITESPARSE_LDL
  // CGInverseOperator + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::LDLOp< DiscreteFunction, LinearOperator >;

    std::string designation(" === LDLOp + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif // HAVE_SUITESPARSE_LDL

#if HAVE_SUITESPARSE_SPQR
  // CGInverseOperator + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::SPQROp< DiscreteFunction, LinearOperator >;

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
  }


  // ISTLInverseOperator< ISTLCGSolver > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::ISTLCGSolver >;

    std::string designation(" === ISTLInverseOperator< ISTLCGSolver > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // ISTLInverseOperator< ISTLBiCGSTABSolver > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::ISTLBiCGSTABSolver >;

    std::string designation(" === ISTLInverseOperator< ISTLBiCGSTABSolver > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // ISTLInverseOperator< ISTLMINRESSolver > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::ISTLMINRESSolver >;

    std::string designation(" === ISTLInverseOperator< ISTLMINRESSolver > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

  // ISTLInverseOperator< ISTLRestartedGMRes > + SparseRowLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLInverseOperator< DiscreteFunction, Dune::Fem::ISTLRestartedGMRes >;

    std::string designation(" === ISTLInverseOperator< ISTLRestartedGMRes > + SparseRowLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }

#if HAVE_SUPERLU
  // ISTLInverseOperator< ISTLSuperLU > + ISTLLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::ISTLLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ISTLSuperLU< DiscreteFunction, LinearOperator >;

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
  }
#endif // HAVE_PETSC

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

#if HAVE_VIENNACL
  // EigenCGInverseOperator + EigenLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::ViennaCLInverseOperator< DiscreteFunction >;

    std::string designation(" === EigenCGInverseOperator + EigenLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif

/*
#if HAVE_AMGXSOLVER
  // EigenBiCGStabInverseOperator + EigenLinearOperator
  {
    using DiscreteFunction  = Dune::Fem::PetscDiscreteFunction< DiscreteSpaceType >;
    using LinearOperator    = Dune::Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction >;
    using InverseOperator   = Dune::Fem::AMGXInverseOperator< DiscreteFunction >;

    std::string designation(" === AMGXInverseOperator + PetscLinearOperator === ");
    pass &= Algorithm< InverseOperator, LinearOperator >::apply( grid, designation, verboseSolver );
  }
#endif
*/

  return pass ? 0 : 1;
}
