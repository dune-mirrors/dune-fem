#ifndef DUNE_FEMPY_GRID_HH
#define DUNE_FEMPY_GRID_HH

#include <cstddef>

#include <memory>
#include <string>

#include <dune/fempy/grid/hierarchical.hh>

namespace Dune
{

  namespace FemPy
  {

    /*!
        @file
        @brief Contains C++ template classes for grids.
        \ingroup Grids
    */

    // LeafGrid
    // --------

    template< class GP >
    struct LeafGrid
    {
      typedef GP GridPart;
      typedef typename GridPart::GridType Grid;

      explicit LeafGrid ( const std::string &dgf )
        : grid_( dgf ),
          gridPart_( new GridPart( *grid_.grid() ) )
      {}

      template< class Mark >
      void adapt ( Mark mark )
      {
        return grid_.adapt( mark );
      }

      void globalRefine ( int level ) { grid_.globalRefine( level ); }

      std::size_t size ( int codim ) const { return gridPart_->indexSet().size( codim ); }

      const std::shared_ptr< Grid > &grid () const { return grid_.grid(); }
      std::shared_ptr< Grid > &grid () { return grid_.grid(); }

      const std::shared_ptr< GridPart > &gridPart () const { return gridPart_; }
      std::shared_ptr< GridPart > &gridPart () { return gridPart_; }

      HierarchicalGrid< Grid > hierarchicalGrid () const { return grid_; }

    private:
      HierarchicalGrid< Grid > grid_;
      std::shared_ptr< GridPart > gridPart_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_HH
