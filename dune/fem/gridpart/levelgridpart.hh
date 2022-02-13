#ifndef DUNE_FEM_GRIDPART_LEVELGRIDPART_HH
#define DUNE_FEM_GRIDPART_LEVELGRIDPART_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/gridview2gridpart.hh>

namespace Dune
{

  namespace Fem
  {

    // LevelGridPart
    // -------------

    template< class Grid >
    class LevelGridPart
      : public GridView2GridPart< typename Grid::LevelGridView, LevelGridPart< Grid > >
    {
      typedef GridView2GridPart< typename Grid::LevelGridView, LevelGridPart< Grid > > BaseType;

    public:
      /** \copydoc Dune::Fem::GridPartInterface::GridType */
      typedef typename BaseType::GridType GridType;

      /** \name Construction
       *  \{
       */

      LevelGridPart ( GridType &grid, int level )
        : BaseType( grid.levelGridView( level ) ),
          grid_( &grid ),
          level_( level )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      using BaseType::grid;

      /** \copydoc Dune::Fem::GridPartInterface::grid */
      GridType &grid () { assert( grid_ ); return *grid_; }

      /** \} */

      /** \brief Return the level that this grid part was created for.
       *  \note This is not an interface method!
       */
      int level() const { return level_; }

    private:
      GridType *grid_;
      int level_;
    };



    namespace GridPartCapabilities
    {

      template< class Grid >
      struct hasGrid< LevelGridPart< Grid > >
      {
        static const bool v = true;
      };

      template< class Grid >
      struct hasSingleGeometryType< LevelGridPart< Grid > >
       : public Dune::Capabilities::hasSingleGeometryType< Grid >
      {};

      template< class Grid >
      struct isCartesian< LevelGridPart< Grid > >
       : public Dune::Capabilities::isCartesian< Grid >
      {};

      template< class Grid, int codim  >
      struct hasEntity< LevelGridPart< Grid >, codim >
       : public Dune::Capabilities::hasEntity< Grid, codim >
      {};

      template< class Grid, int codim  >
      struct canCommunicate< LevelGridPart< Grid >, codim >
       : public Dune::Capabilities::canCommunicate< Grid, codim >
      {};

      template< class Grid >
      struct isConforming< LevelGridPart< Grid > >
      {
        static const bool v = Dune::Capabilities::isLevelwiseConforming< Grid >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_LEVELGRIDPART_HH
