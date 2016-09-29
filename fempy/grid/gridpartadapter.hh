#ifndef DUNE_FEMPY_GRID_GRIDPARTADAPTER_HH
#define DUNE_FEMPY_GRID_GRIDPARTADAPTER_HH

#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/common/gridview2gridpart.hh>

namespace Dune
{

  namespace FemPy
  {

    // GridPartAdapter
    // ---------------

    template< class GridView >
    class GridPartAdapter
      : public Fem::GridView2GridPart< GridView, GridPartAdapter< GridView > >
    {
      typedef GridPartAdapter< GridView > This;
      typedef Fem::GridView2GridPart< GridView, GridPartAdapter< GridView > > Base;

    public:
      typedef typename Base::GridType GridType;

      explicit GridPartAdapter ( const GridView &gridView ) : Base( gridView ) {}

      const GridType &grid () const { return static_cast< GridView >( *this ).grid(); }
      GridType &grid () { return const_cast< GridType & >( static_cast< GridView >( *this ).grid() ); }

      int level () const { DUNE_THROW( NotImplemented, "GridPartAdapter cannot provide level information" ); }
    };

  } // namespace FemPy



  namespace Fem
  {

    namespace GridPartCapabilities
    {

      template< class GridView >
      struct hasGrid< FemPy::GridPartAdapter< GridView > >
      {
        static const bool v = true;
      };

      template< class GridView >
      struct hasSingleGeometryType< FemPy::GridPartAdapter< GridView > >
        : public Dune::Capabilities::hasSingleGeometryType< typename GridView::Grid >
      {};

      template< class GridView >
      struct isCartesian< FemPy::GridPartAdapter< GridView > >
        : public Dune::Capabilities::isCartesian< typename GridView::Grid >
      {};

      template< class GridView, int codim >
      struct hasEntity< FemPy::GridPartAdapter< GridView >, codim >
        : public Dune::Capabilities::hasEntity< typename GridView::Grid, codim >
      {};

      template< class GridView, int codim >
      struct canCommunicate< FemPy::GridPartAdapter< GridView >, codim >
        : public Dune::Capabilities::canCommunicate< typename GridView::Grid, codim >
      {};

      template< class GridView >
      struct isConforming< FemPy::GridPartAdapter< GridView > >
      {
        static const bool v = GridView::conforming;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_GRIDPARTADAPTER_HH
