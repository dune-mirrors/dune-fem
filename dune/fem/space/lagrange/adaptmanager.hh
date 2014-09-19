#ifndef DUNE_FEM_SPACE_LAGRANGE_ADAPTMANAGER_HH
#define DUNE_FEM_SPACE_LAGRANGE_ADAPTMANAGER_HH

// dune-fem includes
#include <dune/fem/space/common/localrestrictprolong.hh>

//local includes
#include "declaration.hh"
#include "restrictprolong.hh"


namespace Dune
{

  namespace Fem
  {
    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class FunctionSpaceType, class GridPartType, int order, template< class > class StorageType >
    class DefaultLocalRestrictProlong< LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, order, StorageType > >
    : public LagrangeLocalRestrictProlong< typename GridPartType::GridType, order >
    {
    public:
      DefaultLocalRestrictProlong ( const LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, order, StorageType > & )
      { }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_ADAPTMANAGER_HH
