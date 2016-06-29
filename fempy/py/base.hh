#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune
{

  namespace FemPy
  {

    // VTKDataType
    // -----------

    enum class VTKDataType { CellData, PointData, CellVector, PointVector };

  }

}
