#ifndef DUNE_FEM_GRIDPART_LEAFGRIDPART_HH
#define DUNE_FEM_GRIDPART_LEAFGRIDPART_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/gridview2gridpart.hh>

namespace Dune
{

  namespace Fem
  {

    // LeafGridPart
    // ------------

    template< class Grid >
    class LeafGridPart
      : public GridView2GridPart< typename Grid::LeafGridView, LeafGridPart< Grid > >
    {
      typedef GridView2GridPart< typename Grid::LeafGridView, LeafGridPart< Grid > > BaseType;

    public:
      /** \copydoc Dune::Fem::GridPartInterface::GridType */
      typedef typename BaseType::GridType GridType;

      /** \name Construction
       *  \{
       */
      explicit LeafGridPart ( GridType &grid )
        : BaseType( grid.leafGridView() ),
          grid_( &grid )
      {}

      LeafGridPart ( const LeafGridPart& other )
        : BaseType( other.grid_->leafGridView() ),
          grid_( other.grid_ )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      using BaseType::grid;

      /** \copydoc Dune::Fem::GridPartInterface::grid */
      GridType &grid () { assert( grid_ ); return *grid_; }

      /** \} */

    private:
      GridType *grid_;
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
