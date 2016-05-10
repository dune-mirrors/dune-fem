#ifndef DUNE_FEMPY_VTK_HH
#define DUNE_FEMPY_VTK_HH

#include <memory>
#include <vector>

#include <dune/fem/io/file/vtkio.hh>

#include <dune/fempy/gridfunction.hh>

namespace Dune
{

  namespace FemPy
  {

    // VTKOutput
    // ---------

    template< class GridPart >
    struct VTKOutput
    {
      explicit VTKOutput ( const LeafGrid< GridPart > &grid )
        : gridPart_( grid.gridPart() ), vtk_( *gridPart_ )
      {}

      template< int R >
      void add ( std::shared_ptr< GridFunction< GridPart, R > > gf )
      {
        gfVector_.push_back( gf );
        vtk_.addCellData( *gf );
      }

      void write ( const char *name ) { vtk_.write( name ); }

    private:
      std::shared_ptr< GridPart > gridPart_;
      Fem::VTKIO< GridPart > vtk_;
      std::vector< std::shared_ptr< GridFunctionBase<GridPart> > > gfVector_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_VTK_HH
