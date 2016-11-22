#include <config.h>

// polynomial order of base functions
const int polOrder = POLORDER;

#include <iostream>
#include <sstream>

#include <dune/fem/gridpart/levelgridpart.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/scalarproducts.hh>


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
  typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;
  typedef Dune::Fem::SlaveDofs< DiscreteFunctionSpaceType, BlockMapperType > SlaveDofsType;
  typedef typename SlaveDofsType :: SingletonKey SlaveDofsKeyType;
  typedef Dune::Fem::SingletonList< SlaveDofsKeyType, SlaveDofsType > SlaveDofsProviderType;
  // test consistency of slavedofs
  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType *space1 = new DiscreteFunctionSpaceType( gridPart );
  DiscreteFunctionSpaceType *space2 = new DiscreteFunctionSpaceType( gridPart );
  SlaveDofsKeyType key1( *space1, space1->blockMapper() );
  SlaveDofsType &slaveDofs1 = SlaveDofsProviderType :: getObject( key1 );
  slaveDofs1.rebuild(*space1);
  std::cout << "number of slave dofs: " << slaveDofs1.size() << std::endl;
  space1->~DiscreteFunctionSpaceType();
  char *tmp = new(space1)char[sizeof(DiscreteFunctionSpaceType)];
  for (int i=0;i<sizeof(DiscreteFunctionSpaceType);++i) tmp[i]=0;
  delete [] tmp;
  SlaveDofsKeyType key2( *space2, space2->blockMapper() );
  SlaveDofsType &slaveDofs2 = SlaveDofsProviderType :: getObject( key2 );
  slaveDofs1.rebuild(*space2);
  std::cout << "number of slave dofs: " << slaveDofs2.size() << std::endl;
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

  std::ostringstream gridName;
  gridName << MyGridType::dimension << "dgrid.dgf";
  Dune::GridPtr< MyGridType > gridptr( gridName.str().c_str() );

  for( int level = 0; level < ml; ++level )
  {
    baseTests( *gridptr, level );
    Dune::Fem::GlobalRefine::apply(*gridptr,1);
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
