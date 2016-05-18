#ifndef DUNE_FEMPY_PY_VTK_HH
#define DUNE_FEMPY_PY_VTK_HH

#include <type_traits>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // VTKDataType
    // -----------

    enum class VTKDataType { CellData, PointData };



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

      pybind11::enum_< VTKDataType > vtkDataType( cls, "DataType" );
      vtkDataType.value( "CellData", VTKDataType::CellData );
      vtkDataType.value( "PointData", VTKDataType::PointData );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_VTK_HH
