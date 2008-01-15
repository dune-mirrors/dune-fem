#ifndef DUNE_LAGRANGEMATRIXSETUP_HH
#define DUNE_LAGRANGEMATRIXSETUP_HH

#include <dune/fem/space/common/gridpartutility.hh>

namespace Dune {

////////////////////////////////////////////////////////////
//
//  Setup of matrix structure 
//
////////////////////////////////////////////////////////////
/**  \brief Setup Matrix structure for Lagrange operators by including
 all lagrange nodes of an element.  
*/
class LagrangeMatrixSetup
{
public:
  //! get number of entries per row for a block matrix, 
  //! i.e. here number of neighboring nodes + 1 
  template <class GridPartImp>
  static inline int stencilSizeEstimate(const GridPartImp& gridPart)        
  {
    return 15;//(GridPartImp :: GridType :: dimension * 2) + 1;
  }

  //! create entries for element and neighbors 
  template <class SpaceImp, 
            class RowMapperType,
            class ColMapperType,
            class MatrixStructureMapImp,
            class DiscreteFunctionType>
  static inline void setup(const SpaceImp& space, 
                           const RowMapperType& rowMapper,
                           const ColMapperType& colMapper,
                           MatrixStructureMapImp& indices,
                           const DiscreteFunctionType*)
  {
    typedef typename SpaceImp :: GridPartType GridPartImp;
    GridPartImp& gridP = const_cast<GridPartImp&> (space.gridPart());
           
    
    typedef typename GridPartNewPartitionType<
      GridPartImp,All_Partition> :: NewGridPartType GridPartType; 

    GridPartType gridPart ( gridP.grid() );

    // define used types 
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridPartType :: template Codim<0> :: IteratorType  IteratorType;

    // clear map 
    indices.clear();

    // only for diagonal 
    IteratorType endit = gridPart.template end<0>(); 
    for(IteratorType it = gridPart.template begin<0>(); it != endit; ++it)
    {
      const EntityType & en = *it;
      // add all column entities to row  
      fill(gridPart,en,rowMapper,colMapper,indices);
    }
  }

protected:
  //! create entries for element and neighbors 
  template <class GridPartImp,
            class EntityImp,
            class RowMapperImp,
            class ColMapperImp,
            class OverlapVectorImp>
  static inline void fill(const GridPartImp& gridPart,
                   const EntityImp& en,
                   const RowMapperImp& rowMapper,
                   const ColMapperImp& colMapper,
                   std::map< int , std::set<int> >& indices)
  {
    // type of local indices storage 
    typedef std::set< int >  LocalIndicesType; 

    const int numBaseFunctions = rowMapper.numDofs(); 

    for(int i=0; i<numBaseFunctions; ++i)
    {
      const int rowIdx = rowMapper.mapToGlobal(en , i );
      LocalIndicesType& localIndices = indices[rowIdx];
      
      // insert diagonal for node  
      localIndices.insert( rowIdx );

      for(int j=0; j<numBaseFunctions; ++j)
      {
        const int colIdx = colMapper.mapToGlobal(en , j );

        // insert column each node  
        localIndices.insert( colIdx );
      }
    }
  }
};

} // end namespace Dune 
#endif
