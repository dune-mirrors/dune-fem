#ifndef DUNE_LAGRANGEMATRIXSETUP_HH
#define DUNE_LAGRANGEMATRIXSETUP_HH

//#include <dune/fem/space/common/gridpartutility.hh>

namespace Dune
{

  /** \class LagrangeMatrixSetup
   *  \brief Setup Matrix structure for Lagrange operators by including all
   *         Lagrange nodes of an element.
   */
  class LagrangeMatrixSetup
  {
  public:
    //! get number of entries per row for a block matrix, 
    //! i.e. here number of neighboring nodes + 1 
    template< class GridPart >
    static inline int stencilSizeEstimate ( const GridPart &gridPart )
    {
      return 15;//(GridPartImp :: GridType :: dimension * 2) + 1;
    }

    //! create entries for element and neighbors 
    template< class Space, class RowMapper, class ColMapper,
              class MatrixStructureMap, class DiscreteFunction >
    static inline void setup ( const Space &space, 
                               const RowMapper &rowMapper,
                               const ColMapper &colMapper,
                               MatrixStructureMap &indices,
                               const DiscreteFunction*)
    {
      typedef typename Space :: IteratorType IteratorType;

      indices.clear();

      const IteratorType end = space.end();
      for( IteratorType it = space.begin(); it != end; ++it )
        fill( space.gridPart(), *it, rowMapper, colMapper, indices );

#if 0
      typedef typename Space :: GridPartType GridPartImp;
      GridPartImp &gridP = const_cast< GridPartImp & >( space.gridPart() );
      
      typedef typename GridPartNewPartitionType< GridPartImp, All_Partition >
        :: NewGridPartType GridPartType;

      GridPartType gridPart ( gridP.grid() );

      // define used types 
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType
        IteratorType;
      typedef typename IteratorType :: Entity EntityType;

      // clear map 
      indices.clear();

      // only for diagonal 
      const IteratorType endit = gridPart.template end< 0 >();
      for( IteratorType it = gridPart.template begin< 0 >(); it != endit; ++it )
      {
        const EntityType &entity = *it;
        // add all column entities to row  
        fill( gridPart, entity, rowMapper, colMapper, indices );
      }
#endif
    }

  protected:
    //! create entries for element and neighbors 
    template< class GridPart, class Entity,
              class RowMapper, class ColMapper >
    static inline void fill ( const GridPart &gridPart,
                              const Entity &entity,
                              const RowMapper &rowMapper,
                              const ColMapper &colMapper,
                              std :: map< int, std :: set< int > > &indices )
    {
      // type of local indices storage 
      typedef std :: set< int >  LocalIndicesType;

      typedef typename RowMapper :: DofMapIteratorType RowDofMapIteratorType;
      typedef typename ColMapper :: DofMapIteratorType ColDofMapIteratorType;

      const RowDofMapIteratorType rowEnd = rowMapper.end( entity );
      RowDofMapIteratorType rowIt = rowMapper.begin( entity );
      for( ; rowIt != rowEnd; ++rowIt )
      {
        LocalIndicesType &localIndices = indices[ rowIt.global() ];

        const ColDofMapIteratorType colEnd = colMapper.end( entity );
        ColDofMapIteratorType colIt = colMapper.begin( entity );
        for( ; colIt != colEnd; ++colIt )
          localIndices.insert( colIt.global() );
      }
#if 0
      const int numBaseFunctions = rowMapper.numDofs( entity ); 

      for( int i=0; i < numBaseFunctions; ++i )
      {
        const int rowIdx = rowMapper.mapToGlobal( entity, i );
        LocalIndicesType &localIndices = indices[ rowIdx ];
        
        // insert diagonal for node  
        localIndices.insert( rowIdx );

        for( int j = 0; j < numBaseFunctions; ++j )
        {
          const int colIdx = colMapper.mapToGlobal( entity, j );

          // insert column each node  
          localIndices.insert( colIdx );
        }
      }
#endif
    }
  };

} // end namespace Dune 
#endif
