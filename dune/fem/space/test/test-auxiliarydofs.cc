#include <config.h>

// polynomial order of base functions
const int polOrder = POLORDER;

#include <iostream>
#include <sstream>

#include <dune/fem/gridpart/levelgridpart.hh>

#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/test/testgrid.hh>


// Type Definitions
// ----------------
typedef Dune::GridSelector::GridType MyGridType;

typedef Dune::Fem::LeafGridPart< MyGridType > GridPartType;

//! type of the function space
typedef Dune::Fem::FunctionSpace< double, double, MyGridType::dimensionworld, 1 > FunctionSpaceType;

//! type of the discrete function space our unkown belongs to
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

// Algorithm to Apply Repeatedly
// -----------------------------

void baseTests ( MyGridType &grid, int level )
{
  typedef typename DiscreteFunctionSpaceType :: AuxiliaryDofsType AuxiliaryDofsType;
  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType *space1 = new DiscreteFunctionSpaceType( gridPart );
  DiscreteFunctionSpaceType *space2 = new DiscreteFunctionSpaceType( gridPart );
  const AuxiliaryDofsType &auxiliaryDofs1 = space1->auxiliaryDofs();
  std::cout << "number of auxiliary dofs: " << auxiliaryDofs1.size() << std::endl;
  space1->~DiscreteFunctionSpaceType();
  char *tmp = new(space1)char[sizeof(DiscreteFunctionSpaceType)];
  for (unsigned int i=0;i<sizeof(DiscreteFunctionSpaceType);++i) tmp[i]=0;
  delete [] tmp;
  const AuxiliaryDofsType &auxiliaryDofs2 = space2->auxiliaryDofs();
  std::cout << "number of auxiliary dofs: " << auxiliaryDofs2.size() << std::endl;
  delete space2;
}


// Main Program
// ------------

int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append parameter
  Dune::Fem::Parameter::append( argc , argv );
  std::string paramFile = "parameter";
  if( argc < 2 )
    std::cerr << "Usage: " << argv[ 0 ] << "<parameter>" << std::endl;
  else
    paramFile = argv[ 1 ];
  Dune::Fem::Parameter::append( paramFile );

  const int ml = Dune::Fem::Parameter::getValue< int >( "lagrangeglobalrefine.maxlevel", 2 );

  std::ostringstream s;
  s << MyGridType::dimension << "dgrid_8.dgf";

  MyGridType &grid = Dune::Fem::TestGrid::grid( s.str() );

  for( int level = 0; level < ml; ++level )
  {
    baseTests( grid, level );
    Dune::Fem::GlobalRefine::apply(grid,1);
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
