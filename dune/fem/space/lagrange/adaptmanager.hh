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

    template< class FunctionSpaceType, class GridPartType, int maxPolOrder, class StorageType >
    class DefaultLocalRestrictProlong< LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, maxPolOrder, StorageType > >
    : public LagrangeLocalRestrictProlong< typename GridPartType::GridType,
                                           // extract maxPolynomialOrder from space because maxPolOrder could be -1
                                           LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, maxPolOrder, StorageType >::maxPolynomialOrder >
    {
    public:
      DefaultLocalRestrictProlong ( const LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, maxPolOrder, StorageType > & )
      { }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_ADAPTMANAGER_HH
