#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>
#include <dune/fem/gridpart/filter/domainfilter.hh>
#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/misc/gridwidth.hh>

#include "failure.hh"
#include "checkseed.hh"
#include "checkgeometry.hh"
#include "checkintersections.hh"
#include "checkgridpart.hh"
#include <dune/fem/test/testgrid.hh>


template< bool UseConsecutiveIndexSet, class HostGridPartType, class FilterType >
void testFilteredGridPart( HostGridPartType& hostGridPart, FilterType& filter )
{
  typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, UseConsecutiveIndexSet > GridPartType;
  GridPartType gridPart( hostGridPart, filter );

  std::cout << "Testing entities" << std::endl;
  testGridPart( gridPart );
  std::cout << std::endl;

  std::cout << "Testing subentities" << std::endl;
  testSubEntities< HostGridPartType::GridType::dimension >( gridPart );
  std::cout << std::endl;

  std::cout << "Testing intersection" << std::endl;
  testIntersectionIterator( gridPart );
  std::cout << std::endl;

  std::cout << "GridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  std::cout << "Testing entity seeds" << std::endl;
  Dune::Fem::CheckEntitySeed< GridPartType >::check( gridPart );
  std::cout << "Testing geometries" << std::endl;
  Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing intersections" << std::endl;
  Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
}

int main( int argc, char ** argv )
{
  Dune::Fem::MPIManager :: initialize( argc, argv );

  try
  {
    // create grid
    auto& grid = Dune::Fem::TestGrid::grid();

    // refine grid
    const int step = Dune::Fem::TestGrid::refineStepsForHalf();
    grid.globalRefine( 2*step );
    grid.loadBalance();

    // create grid part
    typedef Dune::GridSelector::GridType GridType;
    typedef Dune::Fem::AdaptiveLeafGridPart< GridType > HostGridPartType;
    HostGridPartType hostGridPart( grid );

    // create radial filter
    typedef Dune::Fem::RadialFilter< GridType::ctype, GridType::dimensionworld > RadialFilterType;
    typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, RadialFilterType > WrapperRadialFilterType;
    WrapperRadialFilterType wrappedRadialFilter( hostGridPart, RadialFilterType::GlobalCoordinateType( 0 ), 0.25 );

    // test FilteredGridPart with radial filter and allowing non consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with radial filter: allow non consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< false, HostGridPartType, WrapperRadialFilterType >( hostGridPart, wrappedRadialFilter );

    // test FilteredGridPart with radial filter and forcing consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with radial filter: force consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< true, HostGridPartType, WrapperRadialFilterType >( hostGridPart, wrappedRadialFilter );

    // create domain filter
    typedef std::vector< int > DomainArrayType;
    DomainArrayType tags( grid.size(0), 0 );
    for( std::size_t i = ( tags.size()/2 ); i < tags.size(); ++i )
      tags[ i ] = 1;
    typedef Dune::Fem::DomainFilter< HostGridPartType, DomainArrayType > DomainFilterType;
    typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, DomainFilterType > WrapperDomainFilterType;
    WrapperDomainFilterType wrapperDomainFilter( hostGridPart, hostGridPart, tags, 1 );

    // test FilteredGridPart with domain filter and allowing non consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with domain filter: allow non consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< false, HostGridPartType, WrapperDomainFilterType >( hostGridPart, wrapperDomainFilter );

  }
  catch(Dune::Exception &e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
