#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_MAPPERSELECTOR_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_MAPPERSELECTOR_HH

#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/combinedspace/combinedmapper.hh>

#include <dune/fem/space/discontinuousgalerkin.hh>

namespace Dune
{

  namespace Fem
  {

    template<class GridPart, class Mapper1, int blockSize1, class Mapper2, int blockSize2 >
    struct CombinedDiscreteFunctionSpaceMapperSelector;

    template<class GridPart, int codim1, int blockSize1, int codim2, int blockSize2 >
    struct CombinedDiscreteFunctionSpaceMapperSelector< GridPart,
      CodimensionMapper< GridPart, codim1 >, blockSize1,
      CodimensionMapper< GridPart, codim2 >, blockSize2 >
    {
      static_assert( codim1 == codim2, "CombinedDiscreteFunctionSpace only implemented for same codim spaces" );

      static const int localBlockSize = blockSize1 + blockSize2;
      typedef CodimensionMapper< GridPart, codim1 > BlockMapperType;

      typedef CodimensionMapperSingletonFactory< GridPart, codim1 > BlockMapperSingletonFactoryType;
      typedef typename BlockMapperSingletonFactoryType::Key BlockMapperKeyType;
      typedef SingletonList< BlockMapperKeyType, BlockMapperType, BlockMapperSingletonFactoryType > BlockMapperProviderType;

      template< class DFunctionSpace1, class DFunctionSpace2 >
      static BlockMapperType* getBlockMapper ( const DFunctionSpace1 &space1, const DFunctionSpace2 &sapce2 )
      {
        return &BlockMapperProviderType::getObject( space1.gridPart() );
      }

      static void deleteBlockMapper ( BlockMapperType *blockMapper )
      {
        BlockMapperProviderType::removeObject( *blockMapper );
      }
    };


    template<class GridPart, class Mapper1, int blockSize1, class Mapper2, int blockSize2 >
    struct CombinedDiscreteFunctionSpaceMapperSelector
    {
      static const int localBlockSize = 1;
      typedef CombinedSpaceMapper< typename GridPart :: GridType, Mapper1, blockSize1, Mapper2, blockSize2 > BlockMapperType;

      template< class DFunctionSpace1, class DFunctionSpace2 >
      static BlockMapperType* getBlockMapper( const DFunctionSpace1 &space1, const DFunctionSpace2 &space2 )
      {
        return new BlockMapperType( space1.gridPart().grid(), space1.blockMapper(), space2.blockMapper() );
      }

      static void deleteBlockMapper( BlockMapperType *blockMapper )
      {
        delete blockMapper;
        blockMapper = nullptr;
      }
    };

  } // namespace Fem

} // namespace Dune
#endif //#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_MAPPERSELECTOR_HH
