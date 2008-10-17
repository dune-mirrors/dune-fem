#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/gridpart.hh>

#include "testgrid.hh"
#include "dfspace.hh"


template< class EntityType, class MapperType >
inline bool checkIterator ( const EntityType &entity,
                            const MapperType &mapper )
{
  typedef typename MapperType :: DofMapIteratorType DofMapIteratorType;

  int count = 0;
  const DofMapIteratorType end = mapper.end( entity );
  for( DofMapIteratorType it = mapper.begin( entity ); it != end; ++it )
  {
    if( it.global() != mapper.mapToGlobal( entity, it.local() ) )
      return false;
    ++count;
  }
  if( count != mapper.numDofs( entity ) )
    return false;

  return true;
}


template< class DiscreteFunctionSpaceType >
inline bool checkMappers ( const DiscreteFunctionSpaceType &space )
{
  typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
  typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;
  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

  typedef typename IteratorType :: Entity EntityType;
  
  enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

  const MapperType &mapper = space.mapper();
  const BlockMapperType &blockMapper = space.blockMapper();

  const IteratorType end = space.end();
  for( IteratorType it = space.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;

    if( !checkIterator( entity, mapper ) )
    {
      std :: cout << "Mapper: DoF Iterator is inconsistent with DoF mapping."
                  << std :: endl;
      return false;
    }
    
    if( !checkIterator( entity, blockMapper ) )
    {
      std :: cout << "BlockMapper: DoF Iterator is inconsistent with DoF mapping."
                  << std :: endl;
      return false;
    }

    const int numDofs = mapper.numDofs( entity );
    const int numDofBlocks = blockMapper.numDofs( entity );
    /*
    std :: cout << "numDofs: " << numDofs
                << ", numDofBlocks: " << numDofBlocks
                << ", blockSize: " << blockSize
                << std :: endl;
    */
    if( numDofs != blockSize * numDofBlocks )
    {
      std :: cerr << "Number of DoFs is inconsistent between DoF mappers."
                  << std :: endl;
      return false;
    }

    for( int block = 0; block < numDofBlocks; ++block )
    {
      const int blockGlobal = blockMapper.mapToGlobal( entity, block );
      for( int i = 0; i < blockSize; ++i )
      {
        const int dof = blockSize * block + i;
        const int dofGlobal = mapper.mapToGlobal( entity, dof );
        /*
        std :: cout << "local DoF: " << dof << ", globalDoF: " << dofGlobal
                    << ", blockMapped: " << (blockSize * blockGlobal + i)
                    << std :: endl;
        */
        if( dofGlobal != blockSize * blockGlobal + i )
        {
          std :: cerr << "DoF mapping inconsistent between DoF mappers."
                      << std :: endl;
          return false;
        }
      }
    }
  }

  return true;
}


using namespace Dune;

typedef HierarchicGridPart< GridType > GridPartType;

typedef TestFunctionSpace FunctionSpaceType;
typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;


int main ( int argc, char *argv[] )
try
{
  MPIManager :: initialize( argc, argv );

  GridType &grid = TestGrid :: grid();
  const int step = TestGrid :: refineStepsForHalf();
  grid.globalRefine( 2*step );

  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  return (checkMappers( discreteFunctionSpace ) ? 0 : 1);
}
catch( Exception e )
{
  std :: cerr << e.what() << std :: endl;
  return 1;
}
