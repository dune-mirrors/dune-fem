#ifndef DUNE_FEM_SOLVER_PARAMETERDOC_CC
#define DUNE_FEM_SOLVER_PARAMETERDOC_CC

#include <config.h>

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <set>
#include <vector>
#include <iomanip>

#include <dune/grid/yaspgrid.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/matrix/colcompspmatrix.hh>
#include <dune/fem/operator/matrix/istlpreconditioner.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#include <dune/fem/solver/petscavailable.hh>
#include <dune/fem/solver/parameter.hh>

namespace Dune
{
  namespace Fem{
    namespace detail {

      template <typename function_t>
      std::set< std::string > addMethods( const std::vector<int>& mth, function_t* fct )
      {
        std::set< std::string > s;
        for( const auto& m : mth )
          s.insert( fct( m ) );
        return s;
      }

      std::pair< std::string, std::string > solverString(const bool havePetsc)
      {
        std::set< std::string > numpy, istl, petsc, all;
        std::set< std::string > preNum, preIstl, prePetsc, preAll;

        typedef Dune::YaspGrid< 2 > GridType;
        typedef Dune::Fem::FunctionSpace< double , double, GridType::dimensionworld, 1 > FunctionSpaceType;
        typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
        typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 > DiscreteFunctionSpaceType;
        typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunction;

        using Dune::Fem::SolverParameter;

        // numpy
        {
          typedef Dune::Fem::KrylovInverseOperator< DiscreteFunction >  InverseOperatorType;
          numpy = addMethods( InverseOperatorType::supportedSolverMethods(), &SolverParameter::solverMethodTable );
          preNum = addMethods( InverseOperatorType::supportedPreconditionMethods(), &SolverParameter::preconditionMethodTable );
          for( const auto& m : numpy )
            all.insert( m );
          for( const auto& m : preNum )
            preAll.insert( m );
        }

#if HAVE_DUNE_ISTL
        // istl
        {
          istl = addMethods( Dune::Fem::ISTLInverseOperatorMethods::supportedSolverMethods(), &SolverParameter::solverMethodTable );
          preIstl = addMethods( Dune::Fem::ISTLPreconditionMethods::supportedPreconditionMethods(), &SolverParameter::preconditionMethodTable );
          for( const auto& m : istl )
            all.insert( m );
          for( const auto& m : preIstl )
            preAll.insert( m );
        }
#endif

        // petsc
        if( havePetsc )
        {
          typedef Dune::Fem::PetscInverseOperatorAvailable InverseOperatorType;
          petsc = addMethods( InverseOperatorType::supportedSolverMethods(), &SolverParameter::solverMethodTable );
          for( const auto& m : petsc )
            all.insert( m );
          prePetsc = addMethods( InverseOperatorType::supportedPreconditionMethods(), &SolverParameter::preconditionMethodTable );
          auto extra = InverseOperatorType::extraPreconditionMethods();
          for( const auto& p : extra )
            prePetsc.insert( p );
          for( const auto& m : prePetsc )
            preAll.insert( m );
        }

        const auto contains = [&](const std::string& m, const std::set< std::string >& s ) -> std::string
        {
          std::string yes(" x ");
          if( m == "ilu" || m == "lu" )
            yes = std::string(" s ");
          auto it = s.find( m );
          if( it == s.end() )
            return std::string("---");
          else
            return yes;
        };

        std::stringstream out;
        out << "------------------------------------------" << std::endl;
        out << "|  Solver  |         Storage             |" << std::endl;
        out << "|   name   |  numpy  |  istl   |  petsc  |" << std::endl;
        out << "|----------|---------|---------|---------|" << std::endl;
        for( const auto & m : all )
        {
          out << "| " << std::setw(8) << std::left << m << " |   " << contains(m,numpy)
                                                        << "   |   " << contains(m,istl)
                                                        << "   |   " << contains(m,petsc) << "   |" << std::endl;
        }
        out << "------------------------------------------" << std::endl;

        std::stringstream pre;
        pre << "-----------------------------------------------" << std::endl;
        pre << "|  Precondition | (x = parallel | s = serial) |" << std::endl;
        pre << "|  method       |  numpy  |   istl  |  petsc  |" << std::endl;
        pre << "|---------------|---------|---------|---------|" << std::endl;
        for( const auto & m : preAll )
        {
          pre << "| " << std::setw(12) << std::left << m << "  |   " << contains(m,preNum)
                                                        << "   |   " << contains(m,preIstl)
                                                        << "   |   " << contains(m,prePetsc) << "   |" << std::endl;
        }
        pre << "-----------------------------------------------" << std::endl;

        return std::make_pair( out.str(), pre.str() );
      }
} // namespace

} // Fem
} // Dune
#endif
