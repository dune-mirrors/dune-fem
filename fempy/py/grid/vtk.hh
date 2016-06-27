#ifndef DUNE_FEMPY_PY_GRID_VTK_HH
#define DUNE_FEMPY_PY_GRID_VTK_HH

#include <type_traits>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/pybind11/pybind11.h>

#include <dune/fempy/py/base.hh>

namespace Dune
{

  namespace FemPy
  {

    // addToVTKWriter
    // --------------

    template< class GridFunction >
    void addToVTKWriter ( const GridFunction &gf, VTKWriter< typename GridFunction::GridPartType::GridViewType > &vtkWriter, VTKDataType dataType )
    {
      VTK::FieldInfo info( gf.name(), VTK::FieldInfo::Type::scalar, GridFunction::RangeType::dimension );
      switch( dataType )
      {
      case VTKDataType::CellData:
        vtkWriter.addCellData( gf, info );
        break;

      case VTKDataType::PointData:
        vtkWriter.addVertexData( gf, info );
        break;

      default:
        DUNE_THROW( InvalidStateException, "Invalid vtk data type" );
      }
    }



    // registerVTKWriter
    // -----------------

    template< class GridView >
    void registerVTKWriter ( pybind11::handle scope, const char *clsName = "VTKWriter" )
    {
      typedef VTKWriter< GridView > Writer;

      pybind11::class_< Writer > cls( scope, clsName );
      cls.def( "write", [] ( Writer &writer, const std::string &name ) { writer.write( name ); } );
      cls.def( "write", [] ( Writer &writer, const std::string &name, int number ) {
            std::stringstream s; s << name << std::setw(5) << std::setfill('0') << number;
            writer.write( s.str() );
          } );

      cls.attr("CellData") = pybind11::cast(VTKDataType::CellData);
      cls.attr("PointData") = pybind11::cast(VTKDataType::PointData);
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_VTK_HH
