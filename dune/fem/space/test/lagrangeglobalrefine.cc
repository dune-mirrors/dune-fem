#include <config.h>

// polynomial order of base functions
const int polOrder = POLORDER;

#include <iostream>
#include <sstream>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/space/lagrangespace/restrictprolong.hh>
#include <dune/fem/space/common/restrictprolongfunction.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/io/parameter.hh>



// Function to Interpolate
// -----------------------

template< class FunctionSpace >
class ExactSolution
: public Dune::Function< FunctionSpace, ExactSolution< FunctionSpace > >
{
  typedef ExactSolution< FunctionSpace > ThisType;
  typedef Dune::Function< FunctionSpace, ThisType > BaseType;

public:
  typedef FunctionSpace FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

public:
  ExactSolution ( FunctionSpaceType &functionSpace )
  : BaseType( functionSpace )
  {}

  void evaluate ( const DomainType &x, RangeType &phi ) const
  {
    phi = 1;
    for( int i = 0; i < DomainType::dimension; ++i )
      // phi[ 0 ] += x[ i ] * x[ i ]; 
      phi[ 0 ] *= sin( M_PI * x[ i ] ); 
  }

  void evaluate ( const DomainType &x,
                  RangeFieldType t,
                  RangeType &phi ) const
  {
    evaluate( x, phi );
  }

  void jacobian( const DomainType &x, JacobianRangeType &Dphi ) const
  {
    Dphi = 1;
    for( int i = 0; i < DomainType::dimension; ++i )
      for( int j = 0; j < DomainType::dimension; ++j )
        // Dphi[ 0 ][ j ] *= ((i != j) ? 1. : 2.*x[i]);
        Dphi[ 0 ][ j ] *= ((i != j) ? sin( M_PI * x[ i ]) : M_PI * cos( M_PI * x[ i ] ));
  }

  void jacobian( const DomainType &x, RangeFieldType t, JacobianRangeType &Dphi ) const
  {
    jacobian( x, Dphi );
  }
};



// Type Definitions
// ----------------
typedef GridSelector::GridType GridType;

typedef Dune::LevelGridPart< GridType > GridPartType;

//! type of the function space
typedef Dune::FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;

//! type of the discrete function space our unkown belongs to
typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

//! type of the discrete function we are using
typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

typedef Dune::DiscreteFunctionAdapter< ExactSolutionType, GridPartType > GridExactSolutionType;



// Algorithm to Apply Repeatedly
// -----------------------------

void algorithm ( GridType &grid, int level )
{
  const unsigned int polOrder = DiscreteFunctionSpaceType::polynomialOrder;
  FunctionSpaceType functionSpace;
  ExactSolutionType fExact( functionSpace );

  GridPartType fatherGrid( grid, level );
  DiscreteFunctionSpaceType fatherSpace( fatherGrid );
  DiscreteFunctionType fatherFunction( "father function", fatherSpace );

  GridExactSolutionType fatherExact( "father level exact solution", fExact, fatherGrid, polOrder+1 );
  Dune::L2Norm< GridPartType > fatherL2norm( fatherGrid );
  Dune::H1Norm< GridPartType > fatherH1norm( fatherGrid );

  Dune::LagrangeInterpolation< DiscreteFunctionType >::interpolateFunction( fatherExact, fatherFunction );
  double fatherL2error = fatherL2norm.distance( fatherExact, fatherFunction );
  double fatherH1error = fatherH1norm.distance( fatherExact, fatherFunction );

  std::cout << "Unknowns on father level: " << fatherSpace.size() << std::endl;
  std::cout << "L2 error on father level: " << fatherL2error << std::endl;
  std::cout << "H1 error on father level: " << fatherH1error << std::endl; 

  GridPartType sonGrid( grid, level+1 );
  DiscreteFunctionSpaceType sonSpace( sonGrid );
  DiscreteFunctionType sonFunction( "son function", sonSpace );

  GridExactSolutionType sonExact( "son level exact solution", fExact, sonGrid, polOrder+1 );
  Dune::L2Norm< GridPartType > sonL2norm( sonGrid );
  Dune::H1Norm< GridPartType > sonH1norm( sonGrid );

  Dune::ProlongFunction< Dune::LagrangeLocalRestrictProlong< GridType, polOrder > > prolongFunction;
  prolongFunction( fatherFunction, sonFunction );
  double sonL2error = sonL2norm.distance( sonExact, sonFunction );
  double sonH1error = sonH1norm.distance( sonExact, sonFunction );

  std::cout << "Unknowns on son level: " << sonSpace.size() << std::endl;
  std::cout << "L2 error on son level: " << sonL2error << std::endl;
  std::cout << "H1 error on son level: " << sonH1error << std::endl; 

  double l2eoc = -log( sonL2error / fatherL2error) / M_LN2;
  double h1eoc = -log( sonH1error / fatherH1error) / M_LN2;

  std::cout << "L2 EOC: " << l2eoc << std :: endl;
  std::cout << "H1 EOC: " << h1eoc << std :: endl;

  std::cout << std::endl;
}



// Main Program
// ------------

int main ( int argc, char **argv )
try
{
  Dune::MPIManager::initialize( argc, argv );

  // append parameter
  Dune::Parameter::append( argc , argv );
  std::string paramFile = "parameter";
  if( argc < 2 )
    std::cerr << "Usage: " << argv[ 0 ] << "<parameter>" << std::endl;
  else 
    paramFile = argv[ 1 ]; 
  Dune::Parameter::append( paramFile );

  const int ml = Dune::Parameter::getValue< int >( "lagrangeglobalrefine.maxlevel", 2 );

  std::ostringstream gridName;
  gridName << dimworld << "dgrid.dgf";
  Dune::GridPtr< GridType > gridptr( gridName.str().c_str() );

  //const int step = DGFGridInfo< GridType >::refineStepsForHalf();
  gridptr->globalRefine( ml );

  for( int level = 0; level < gridptr->maxLevel(); ++level )
    algorithm( *gridptr, level );

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
