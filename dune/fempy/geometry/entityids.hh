#ifndef DUNE_FEMPY_GEOMETRY_ENTITYIDS_HH
#define DUNE_FEMPY_GEOMETRY_ENTITYIDS_HH

#include <dune/fem/misc/boundaryidprovider.hh>

namespace Dune
{
  namespace FemPy
  {

    template< class GridPart, class Intersection>
    inline static int boundaryId ( const Intersection &intersection )
    {
      return Dune::Fem::BoundaryIdProvider< typename GridPart::GridType > ::
             boundaryId( intersection );
    }

  }
}

#endif // DUNE_FEMPY_GEOMETRY_ENTITYIDS_HH
