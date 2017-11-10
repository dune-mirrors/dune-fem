#include <config.h>

#include <iostream>
#include <cmath>

#include "../../test/testgrid.hh"
#include "../../test/exactsolution.hh"

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/common/interpolate.hh>

int main ( int argc, char **argv )
{
  //initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // create grid
  typedef Dune::GridSelector::GridType GridType;
  GridType &grid = Dune::Fem::TestGrid::grid();

  // create grid part
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );

  // function space type
  typedef typename GridPartType::ctype DomainFieldType;
  constexpr int dimDomain = GridPartType::dimensionworld;
  constexpr int dimRange = DIMRANGE;
  constexpr int polOrder = POLORDER;
  typedef Dune::Fem::FunctionSpace< DomainFieldType, double, dimDomain, dimRange > FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  // create discrete function space
  typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

  // create discrete functions
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  DiscreteFunctionType dfLocalAdapted( "df local adapted", discreteFunctionSpace );
  dfLocalAdapted.clear();
  DiscreteFunctionType dfGridAdapted( "df grid adapted", discreteFunctionSpace );
  dfGridAdapted.clear();

  // create local function adapter
  typedef Dune::Fem::LocalAnalyticalFunctionBinder<DiscreteFunctionSpaceType> LocalAnalyticalFunctionType;
  LocalAnalyticalFunctionType localAnalyticalFunction(
      [&](const DomainType &x,double t,const EntityType& entity)
      {
        RangeType ret(1);
        for(int r = 0; r < RangeType :: dimension; ++r )
          for( int i = 0; i < DomainType :: dimension; ++i )
            ret[ r ] += pow(sin( M_PI * x[ i ] ),double(r+1));
        #if defined USE_COMPLEX
        ret*= std::complex<double>( 1 , -2. );
        #endif
        return ret;
      });
  typedef Dune::Fem::LocalFunctionAdapter<LocalAnalyticalFunctionType> LocalAdaptedFunctionType;
  LocalAdaptedFunctionType localAdapted("local adapted function",std::move(localAnalyticalFunction),gridPart,5);

  // interpolate local adpated function over discrete function
  Dune::Fem::interpolate(localAdapted,dfLocalAdapted);

  // create grid function adapter
  typedef Dune::Fem::ExactSolution< FunctionSpaceType > GridAnalyticalFunctionType;
  GridAnalyticalFunctionType gridAnalyticalFunction;
  typedef Dune::Fem::GridFunctionAdapter< GridAnalyticalFunctionType, GridPartType > GridAdaptedFunctionType;
  GridAdaptedFunctionType gridAdapted( "grid adapted function", gridAnalyticalFunction, gridPart, 5 );

  // interpolate grid adpated function over discrete function
  Dune::Fem::interpolate( gridAdapted, dfGridAdapted );

  // compare the 2 adapted functions
  bool areDifferent(false);
  auto it(dfGridAdapted.dbegin());
  for(const auto& dof:dofs(dfLocalAdapted))
  {
    if(std::abs(dof-*it)>1.e-12)
      areDifferent=true;
    ++it;
  }
  if(areDifferent)
    std::cout<<"ERROR: the 2 adapted functions are different!"<<std::endl;
  else
    std::cout<<"The 2 adapted functions are equal!"<<std::endl;

  return 0;
}
