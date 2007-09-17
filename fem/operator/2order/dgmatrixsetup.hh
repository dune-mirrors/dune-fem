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
  //! create entries for element and neighbors 
  template <class RowSpaceType,
            class ColSpaceType,
            class MatrixStructureMapImp,
            class OverlapVectorImp>
  static inline void setup(const RowSpaceType& rowSpace,
                           const ColSpaceType& colSpace,
                           MatrixStructureMapImp& indices,
                           OverlapVectorImp& overlapRows)
  {
    typedef typename GridPartNewPartitionType<typename RowSpaceType ::
      GridPartType,All_Partition> :: NewGridPartType GridPartType; 

    GridPartType gridPart ( const_cast<RowSpaceType&>
        (rowSpace).gridPart().grid());

    // define used types 
    typedef typename RowSpaceType :: GridType GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridPartType :: template Codim<0> :: IteratorType  IteratorType;

    typedef typename RowSpaceType :: IndexSetType RowIndexSetType;
    const RowIndexSetType& rowSet = rowSpace.indexSet();

    typedef typename ColSpaceType :: IndexSetType ColIndexSetType;
    const ColIndexSetType& colSet = colSpace.indexSet();

    // clear map 
    indices.clear();

    // we need All_Partition here to insert overlap entities 
    // only for diagonal 
    IteratorType endit = gridPart.template end<0>(); 
    for(IteratorType it = gridPart.template begin<0>(); it != endit; ++it)
    {
      const EntityType & en = *it;
      // add all column entities to row  
      fill(gridPart,en,rowSet,colSet,indices,overlapRows);
    }
  }

protected:
  //! create entries for element and neighbors 
  template <class GridPartImp,
            class EntityImp,
            class RowSetImp,
            class ColSetImp,
            class OverlapVectorImp>
  static inline void fill(const GridPartImp& gridPart,
                   const EntityImp& en,
                   const RowSetImp& rowSet,
                   const ColSetImp& colSet,
                   std::map< int , std::set<int> >& indices,
                   OverlapVectorImp& overlapRows)
  {
    // get index for entity 
    const int elRowIndex = rowSet.index(en);

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
        const int nbColIndex = colSet.index( nb );
        const int nbRowIndex = rowSet.index( nb );

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
          const int elColIndex = colSet.index( en );
          nbIndices.insert( nbColIndex );
          nbIndices.insert( elColIndex );  
        }
      }
    }
  }
};

} // end namespace Dune 
#endif
