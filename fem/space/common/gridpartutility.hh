#ifndef DUNE_GRIDPARTUTILITY_HH
#define DUNE_GRIDPARTUTILITY_HH

//- System includes

//- Dune includes
#include <dune/grid/common/gridpart.hh>

/** @file
  @author Robert Kloefkorn
  @brief Provides other partition types of the same grid part 
*/
/*! @addtogroup GridPart
    @ingroup Grid

    Grid Part Utilities.
*/

namespace Dune {
  /** 
   * @addtogroup GridPart
   *
   * @{ 
   */

  //! \brief Helper class to get GridPart with different partition type
  template <class GridPartImp, PartitionIteratorType piType> 
  struct GridPartNewPartitionType;

  //! \brief Helper class to get GridPart with different partition type
  template <class GridImp, PartitionIteratorType old_pitype, 
            template <class, PartitionIteratorType> class GridPartImp, 
            PartitionIteratorType pitype> 
  struct GridPartNewPartitionType<GridPartImp<GridImp, old_pitype> , pitype>
  {
    //! type of grid part with new partition type 
    typedef GridPartImp<GridImp,pitype> NewGridPartType;
  };
  /** @} */

} // end namespace Dune
#endif
