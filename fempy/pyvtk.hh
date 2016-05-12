#ifndef DUNE_FEMPY_PYVTK_HH
#define DUNE_FEMPY_PYVTK_HH

#include <type_traits>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // PyVTKWriteable
    // --------------

    template< class GridView >
    struct PyVTKWriteable
    {
      typedef VTKWriter< GridView > Writer;

    private:
      typedef VTK::FieldInfo Info;

      struct Interface
      {
        virtual ~Interface () = default;
        virtual void addAsCellData ( Writer &writer ) = 0;
        virtual void addAsPointData ( Writer &writer ) = 0;
      };

      template< class View >
      struct Implementation
      {
        Implementation ( View &&view, Info info ) : view_( std::move( view ) ), info_( std::move( info ) ) {}

        virtual void addAsCellData ( Writer &writer ) { writer.addCellData( view_, info_ ); }
        virtual void addAsPointData ( Writer &writer ) { writer.addVertexData( view_, info_ ); }

      private:
        View view_;
        Info info_;
      };

    public:
      template< class GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value, int > = 0 >
      PyVTKWriteable ( const GF &gf )
      {
        typedef GridFunctionView< GF > View;
        Info info( gf.name(), Info::Type::scalar, GF::RangeType::dimension );
        impl_.reset( new Implementation< View >( gf, info ) );
      }

      void addAsCellData ( Writer &writer ) { impl_->addAsCellData( writer ); }
      void addAsPointData ( Writer &writer ) { impl_->addAsPointData( writer ); }

    private:
      std::unique_ptr< Interface > impl_;
    };



    // registerVTKWriter
    // -----------------

    template< class GridView >
    void registerVTKWriter ( pybind11::handle scope, const char *clsName = "VTKWriter" )
    {
      typedef VTKWriter< GridView > Writer;
      typedef PyVTKWriteable< GridView > Writeable;

      pybind11::class_< Writeable > clsWriteable( scope, "VTKWriteable" );

      pybind11::class_< Writer > cls( scope, clsName );
      cls.def( "addCellData", [] ( Writer &writer, Writeable &writeable ) { writeable.addAsCellData( writer ); } );
      cls.def( "addPointData", [] ( Writer &writer, Writeable &writeable ) { writeable.addAsPointData( writer ); } );
      cls.def( "write", [] ( Writer &writer, const std::string &name ) { writer.write( name ); } );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PYVTK_HH
