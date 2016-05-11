#ifndef DUNE_FEMPY_PYVTK_HH
#define DUNE_FEMPY_PYVTK_HH

#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    template< class VTKWriter >
    void registerVTKWriter ( pybind11::handle scope, const char *clsName = "VTKWriter" )
    {
      pybind11::class_< VTKWriter > cls( scope, clsName );
      cls.def( "write", [] ( VTKWriter &vtk, const std::string &name ) { vtk.write( name ); } );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PYVTK_HH
