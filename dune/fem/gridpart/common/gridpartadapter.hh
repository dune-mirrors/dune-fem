#ifndef DUNE_FEM_GRID_GRIDPARTADAPTER_HH
#define DUNE_FEM_GRID_GRIDPARTADAPTER_HH

#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/common/gridview2gridpart.hh>

namespace Dune
{

  namespace Fem
  {

    // GridPartAdapter
    // ---------------

    template< class GridView >
    class GridPartAdapter
      : public Fem::GridView2GridPart< GridView, GridPartAdapter< GridView >, false >
    {
      typedef GridPartAdapter< GridView > This;
      typedef Fem::GridView2GridPart< GridView, GridPartAdapter< GridView >, false > Base;

    public:
      typedef GridView GridViewType;
      typedef typename Base::GridType GridType;

      explicit GridPartAdapter ( const GridView &gridView ) : Base( gridView ) {}
    };


    namespace GridPartCapabilities
    {

      template< class GridView >
      struct hasGrid< GridPartAdapter< GridView > >
      {
        static const bool v = true;
      };

      template< class GridView >
      struct hasSingleGeometryType< GridPartAdapter< GridView > >
        : public Dune::Capabilities::hasSingleGeometryType< typename GridView::Grid >
      {};

      template< class GridView >
      struct isCartesian< GridPartAdapter< GridView > >
        : public Dune::Capabilities::isCartesian< typename GridView::Grid >
      {};

      template< class GridView, int codim >
      struct hasEntity< GridPartAdapter< GridView >, codim >
        : public Dune::Capabilities::hasEntity< typename GridView::Grid, codim >
      {};

      template< class GridView, int codim >
      struct canCommunicate< GridPartAdapter< GridView >, codim >
        : public Dune::Capabilities::canCommunicate< typename GridView::Grid, codim >
      {};

      template< class GridView >
      struct isConforming< GridPartAdapter< GridView > >
      {
        static const bool v = GridView::conforming;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRID_GRIDPARTADAPTER_HH
