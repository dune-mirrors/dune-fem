#include <config.h>

#include <cassert>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/geogridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/space/lagrange.hh>

#include "./failure.hh"
#include "./checkseed.hh"
#include "./checkgeometry.hh"
#include "./checkindexset.hh"
#include "./checkintersections.hh"
#include "../../test/testgrid.hh"

template< class GridPartType >
void testGridPart( const GridPartType & gridPart )
{
  typedef typename GridPartType::IndexSetType IndexSetType;
  typedef typename GridPartType::IndexSetType::IndexType IndexType;
  const IndexSetType & indexSet = gridPart.indexSet();
  IndexType maxIndex = 0;
  int count = 0;
  std::vector< bool > isConsecutiveIndex(GridPartType::dimension+1,true);
  std::vector< std::vector< bool > > consecutiveIndex(GridPartType::dimension+1);
  for (int c=0;c<=GridPartType::dimension;++c)
    consecutiveIndex[c].resize( indexSet.size(c), false );

  typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< 0 >();
  for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
  {
    ++count;
    IndexType index = indexSet.index(*it);
    maxIndex = std::max( index, maxIndex);

    if (index >= consecutiveIndex[0].size()) isConsecutiveIndex[0] = false;
    else consecutiveIndex[0][index] = true;
    for (int c=0;c<=GridPartType::dimension;++c)
    {
      int nSubEn = Dune::ReferenceElements< typename GridPartType::GridType::ctype, GridPartType::dimension >::
          general( it->type() ).size(c);
      for (int i=0;i<nSubEn;++i)
      {
        IndexType index = indexSet.subIndex(*it,i,c);
        if (index >= consecutiveIndex[c].size()) isConsecutiveIndex[c] = false;
        else consecutiveIndex[c][index] = true;
      }
    }
  }

  std::cout << "entities visited: " << count << std::endl;
  
  std::cout << "entities in index set: " << indexSet.size( 0 ) << std::endl;
  std::cout << "maximum value in index set: " << maxIndex << std::endl;

  for (int c=0;c<=GridPartType::dimension;++c)
  {
    std::cout << "index set for codim " << c; 
    if ( !isConsecutiveIndex[c] ) 
    {
      std::cout << " is not consecutive: too large index encountered" << std::endl;
      continue;
    }
    for (unsigned int i=0;i<consecutiveIndex[c].size();++i)
    {
      if (!consecutiveIndex[c][i]) 
      {
        isConsecutiveIndex[c] = false;
        break;
      }
    }
    if ( !isConsecutiveIndex[c] ) 
    {
      std::cout << " is not consecutive: hole encountered" << std::endl;
    }
    else
    {
      std::cout << " is consecutive" << std::endl;
    }
  }
}

template< int codim, class GridPartType >
void testSubEntities( const GridPartType & gridPart )
{
  typedef typename GridPartType::IndexSetType IndexSetType;
  typedef typename GridPartType::IndexSetType::IndexType IndexType;
  const IndexSetType & indexSet = gridPart.indexSet();
  IndexType maxIndex = 0;
  int count = 0;
  typedef typename GridPartType::template Codim< codim >::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< codim >();
  for( IteratorType it = gridPart.template begin< codim >(); it != end; ++it )
  {
    ++count;
    IndexType index = indexSet.index(*it);
    maxIndex = std::max( index, maxIndex);
  }

  std::cout << "codim " << codim << " subentities visited: " << count << std::endl;

  std::cout << "entities in index set: " << indexSet.size( codim ) << std::endl;
  std::cout << "maximum value in index set: " << maxIndex << std::endl;
}

template< class GridPartType >
void testIntersectionIterator( const GridPartType & gridPart )
{
  std::vector<int> index( gridPart.indexSet().size(0), 0 );
  typedef typename GridPartType::template Codim< 0 >::template Partition<Dune::All_Partition>::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< 0,Dune::All_Partition >();
  for( IteratorType it = gridPart.template begin< 0,Dune::All_Partition >(); it != end; ++it )
    index[ gridPart.indexSet().index( * it ) ] = 1;
  for( IteratorType it = gridPart.template begin< 0,Dune::All_Partition >(); it != end; ++it )
  {
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    const IntersectionIteratorType iend = gridPart.iend( *it );
    for ( IntersectionIteratorType inter = gridPart.ibegin( *it );
          inter != iend; ++inter )
    {
      if (inter->neighbor())
      {
        typename GridPartType::IndexSetType::IndexType nbIndex = gridPart.indexSet().index( *(inter->outside()) );
        if ( nbIndex >= index.size() )
        {
          std::cout << "An index on neighbor is too large" << std::endl;
          continue;
        }
        if ( index[ nbIndex ] == 0 )
        {
          std::cout << "A neighbor is not part of the gridPart" << std::endl;
          continue;
        }
      }
    }
  }
}


template< class FunctionSpace >
class Identity
: public Dune::Fem::Function< FunctionSpace, Identity< FunctionSpace > >
{
  typedef Dune::Fem::Function< FunctionSpace, Identity< FunctionSpace > > BaseType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = x;
  }
};


template< class GridPart, class LocalFunction >
void testExchangeGeometry ( const GridPart &gridPart, LocalFunction &localFunction )
{
  typedef typename GridPart::template Codim< 0 >::IteratorType IteratorType;
  typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPart::template Codim< 0 >::EntityPointerType EntityPointerType;
  typedef typename GridPart::template Codim< 0 >::GeometryType GeometryType;

  typedef Dune::ReferenceElement< typename GridPart::ctype, GridPart::dimension > RefElementType;
  typedef Dune::ReferenceElements< typename GridPart::ctype, GridPart::dimension > RefElementsType;


  const IteratorType end = gridPart.template end< 0 >();
  for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
  {
    const EntityType &entity = *it;

    const RefElementType &refElement = RefElementsType::general( entity.type() );
    localFunction.init( entity.impl().hostEntity() );
    for( int i = 0; i < refElement.size( GridPart::dimension ); ++i )
    {
      for( int k = 0; k < GridPart::dimensionworld; ++k )
        localFunction[ i*GridPart::dimensionworld + k ] = refElement.position( i, GridPart::dimension )[ k ];
    }

    const EntityPointerType xchgEntityPointer = gridPart.exchangeGeometry( entity, localFunction );
    const EntityType &xchgEntity = *xchgEntityPointer;

    const GeometryType &xchgGeometry = xchgEntity.geometry();
    if( xchgGeometry.type() != entity.type() )
      DUNE_THROW( Dune::InvalidStateException, "exchangeGeometry returns wrong geometry type." );
    if( (xchgGeometry.center() - refElement.position( 0, 0 )).two_norm() > 1e-8 )
    {
      std::cerr << "exchangeGeometry returns wrong center: " << xchgGeometry.center() << std::endl;
      std::cerr << "(real geometry center is: " << entity.geometry().center() << ")." << std::endl;
    }
  }
};


typedef Dune::GridSelector::GridType GridType;
//typedef Dune::Fem::AdaptiveLeafGridPart< GridType > HostGridPartType;
typedef Dune::Fem::LeafGridPart< GridType > HostGridPartType;
typedef Dune::Fem::FunctionSpace< GridType::ctype, GridType::ctype, GridType::dimensionworld, GridType::dimensionworld > CoordFunctionSpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< CoordFunctionSpaceType, HostGridPartType, 1 > DiscreteCoordFunctionSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteCoordFunctionSpaceType > CoordFunctionType;
typedef Dune::Fem::GeoGridPart< CoordFunctionType > GridPartType;


int main ( int argc, char ** argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // create grid
  GridType &grid = Dune::TestGrid::grid();

  // refine grid
  const int step = Dune::TestGrid::refineStepsForHalf();
  grid.globalRefine( 2*step );
  grid.loadBalance();

  // create grid part
  HostGridPartType hostGridPart( grid );
  DiscreteCoordFunctionSpaceType coordFunctionSpace( hostGridPart );
  CoordFunctionType coordFunction( "coordinate function", coordFunctionSpace );
  typedef Identity< CoordFunctionSpaceType > IdentityType;
  Dune::Fem::LagrangeInterpolation< IdentityType, CoordFunctionType >::interpolateFunction( IdentityType(), coordFunction );
  GridPartType gridPart( coordFunction );

  // run test
  std::cout << "test using codim=0 iterator" << std::endl;
  testGridPart( gridPart );
  std::cout << std::endl;
  std::cout << "test using codim=dimension iterator" << std::endl;
  testSubEntities< GridType::dimension >( gridPart );

  std::cout << std::endl;
  std::cout << "gridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

  Dune::Fem::TemporaryLocalFunction< DiscreteCoordFunctionSpaceType > tlf( coordFunctionSpace );
  testExchangeGeometry( gridPart, tlf );

  std::cout << std::endl << std::endl;

  // check entity seed
  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  Dune::Fem::CheckEntitySeed< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  Dune::Fem::CheckIndexSet< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
