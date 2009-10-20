#ifndef DUNE_FEM_ENTITYCOMMHELPER_HH
#define DUNE_FEM_ENTITYCOMMHELPER_HH

#include <dune/grid/common/grid.hh>

namespace Dune
{

  template< InterfaceType interface >
  struct EntityCommHelper;


  template<>
  struct EntityCommHelper< InteriorBorder_InteriorBorder_Interface >
  {
    inline static bool send ( const PartitionType p )
    {
      return (p == BorderEntity);
    }

    inline static bool receive ( const PartitionType p )
    {
      return (p == BorderEntity);
    }
  };

  
  template<>
  struct EntityCommHelper< InteriorBorder_All_Interface >
  {
    inline static bool send ( const PartitionType p )
    {
      return (p == InteriorEntity) || (p == BorderEntity);
    }

    inline static bool receive ( const PartitionType p )
    {
      return (p != InteriorEntity);
    }
  };


  template<>
  struct EntityCommHelper< All_All_Interface >
  {
    inline static bool send ( const PartitionType p )
    {
      return true;
    }

    inline static bool receive ( const PartitionType p )
    {
      return true;
    }
  };

}

#endif
