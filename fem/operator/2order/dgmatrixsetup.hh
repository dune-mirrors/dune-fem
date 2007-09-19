#ifndef DUNE_DGMATRIXSETUP_HH
#define DUNE_DGMATRIXSETUP_HH

#include <dune/fem/space/common/gridpartutility.hh>

namespace Dune {

////////////////////////////////////////////////////////////
//
//  Setup of matrix structure 
//
////////////////////////////////////////////////////////////
/**  \brief Setup Matrix structure for DG operators by including
 * elements and it's neighbors. 
*/
class ElementAndNeighbors
{
public:
  //! get number of entries per row for a block matrix, 
  //! i.e. here number of neighbors + 1
  template <class GridPartImp> 
  static inline int stencilSizeEstimate(const GridPartImp& gridPart) 
  {
    return (GridPartImp :: GridType :: dimension * 2) + 1; 
  }

  //! create entries for element and neighbors 
  template <class GridPartImp,    
            class RowMapperType,
            class ColMapperType,
            class MatrixStructureMapImp,
            class OverlapVectorImp>
  static inline void setup(const GridPartImp& gP,    
                           const RowMapperType& rowMapper,
                           const ColMapperType& colMapper,
                           MatrixStructureMapImp& indices,
                           OverlapVectorImp& overlapRows)
  {
    typedef typename GridPartNewPartitionType<
      GridPartImp,All_Partition> :: NewGridPartType GridPartType;    

    GridPartType gridPart ( const_cast<GridPartImp&> (gP).grid());

    // define used types 
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridPartType :: template Codim<0> :: IteratorType  IteratorType;

    // clear map 
    indices.clear();

    // we need All_Partition here to insert overlap entities 
    // only for diagonal 
    IteratorType endit = gridPart.template end<0>(); 
    for(IteratorType it = gridPart.template begin<0>(); it != endit; ++it)
    {
      const EntityType & en = *it;
      // add all column entities to row  
      fill(gridPart,en,rowMapper,colMapper,indices,overlapRows);
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
                   std::map< int , std::set<int> >& indices,
                   OverlapVectorImp& overlapRows)
  {
    assert( rowMapper.numDofs () == 1 );
    // get index for entity 
    const int elRowIndex = rowMapper.mapToGlobal( en, 0 ); 

    // type of local indices storage 
    typedef std::set< int >  LocalIndicesType; 
    LocalIndicesType& localIndices = indices[elRowIndex];

    // insert diagonal for each element 
    localIndices.insert( elRowIndex );

    // if entity is not interior, insert into overlap entries 
    if(en.partitionType() != InteriorEntity)
    {
      overlapRows.push_back( elRowIndex );
      return ;
    }

    // insert neighbors 
    typedef typename GridPartImp:: GridType :: template Codim<0>::EntityPointer EntityPointerType; 
    typedef typename GridPartImp:: IntersectionIteratorType IntersectionIteratorType;
    IntersectionIteratorType endnit = gridPart.iend(en);
    for(IntersectionIteratorType nit = gridPart.ibegin(en);
        nit != endnit; ++nit)
    {
      if(nit.neighbor())
      {
        // get neighbor 
        EntityPointerType ep = nit.outside();
        const EntityImp& nb = *ep;

        // get index of neighbor 
        const int nbColIndex = colMapper.mapToGlobal( nb , 0 );
        const int nbRowIndex = rowMapper.mapToGlobal( nb , 0 );

        // check whether to insert now 
        bool insertHere = (elRowIndex < nbRowIndex);
#if HAVE_MPI 
        // check partition type 
        if( nb.partitionType() != InteriorEntity )
        {
          insertHere = true;
          overlapRows.push_back( nbRowIndex );
        }
#endif
        // insert pair 
        if( insertHere )
        {
          // insert neihgbor 
          localIndices.insert( nbColIndex );

          // insert symetric part with swaped row-col
          LocalIndicesType& nbIndices = indices[nbRowIndex];
          const int elColIndex = colMapper.mapToGlobal( en , 0 );
          nbIndices.insert( nbColIndex );
          nbIndices.insert( elColIndex );  
        }
      }
    }
  }
};

} // end namespace Dune 
#endif
