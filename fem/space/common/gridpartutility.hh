#ifndef DUNE_GRIDPARTUTILITY_HH
#define DUNE_GRIDPARTUTILITY_HH

//- System includes

//- Dune includes
#include <dune/grid/common/gridpart.hh>

/** @file
  @brief Provides other partition types of the same grid part 
*/

namespace Dune {
  /** 
   * @addtogroup GridPart
   *
   * @{ 
   */

  //! @ingroup GridPart
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
