#ifndef DUNE_FEM_GRIDPART_LEAFGRIDPART_HH
#define DUNE_FEM_GRIDPART_LEAFGRIDPART_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/gridview2gridpart.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>

namespace Dune
{

  namespace Fem
  {

    // LeafGridPart
    // ------------

    template< class Grid >
    class LeafGridPart
      : public GridView2GridPart< typename Grid::LeafGridView, WrappedLeafIndexSet< Grid >, LeafGridPart< Grid > >
    {
      typedef GridView2GridPart< typename Grid::LeafGridView, WrappedLeafIndexSet< Grid >, LeafGridPart< Grid > > BaseType;

    public:
      /** \copydoc Dune::Fem::GridPartInterface::GridType */
      typedef typename BaseType::GridType GridType;

      /** \copydoc Dune::Fem::GridPartInterface::GridType */
      typedef typename BaseType::IndexSetType IndexSetType;

      /** \name Construction
       *  \{
       */

      explicit LeafGridPart ( GridType &grid )
        : BaseType( grid.leafGridView() ),
          grid_( grid ),
          indexSet_( grid )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::GridPartInterface::grid */
      GridType &grid () { return grid_; }

      /** \copydoc Dune::Fem::GridPartInterface::grid */
      const GridType &grid () const { return grid_; }

      /** \copydoc Dune::Fem::GridPartInterface::indexSet */
      const IndexSetType &indexSet () const { return indexSet_; }

      /** \copydoc Dune::Fem::GridPartInterface::level */
      int level () const { return grid().maxLevel(); }

      /** \} */

    private:
      GridType &grid_;
      IndexSetType indexSet_;
    };



    namespace GridPartCapabilities
    {

      template< class Grid >
      struct hasGrid< LeafGridPart< Grid > >
      {
        static const bool v = true;
      };

      template< class Grid >
      struct hasSingleGeometryType< LeafGridPart< Grid > >
       : public Dune::Capabilities::hasSingleGeometryType< Grid >
      {};

      template< class Grid >
      struct isCartesian< LeafGridPart< Grid > >
       : public Dune::Capabilities::isCartesian< Grid >
      {};

      template< class Grid, int codim  >
      struct hasEntity< LeafGridPart< Grid >, codim >
       : public Dune::Capabilities::hasEntity< Grid, codim >
      {};

      template< class Grid >
      struct isParallel< LeafGridPart< Grid > >
       : public Dune::Capabilities::isParallel< Grid >
      {};

      template< class Grid, int codim  >
      struct canCommunicate< LeafGridPart< Grid >, codim >
       : public Dune::Capabilities::canCommunicate< Grid, codim >
      {};

      template< class Grid >
      struct isConforming< LeafGridPart< Grid > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< Grid >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_LEAFGRIDPART_HH
