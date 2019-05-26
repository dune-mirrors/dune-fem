#ifndef DUNE_FEM_GRIDPART_COMMON_UNIMPLEMENTEDGEOMETRY_HH
#define DUNE_FEM_GRIDPART_COMMON_UNIMPLEMENTEDGEOMETRY_HH

#include <string>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  // UnimplementedGeometry
  // ---------------------

  template< class ct, int mydim, int cdim >
  struct UnimplementedGeometry
  {
    static const int mydimension = mydim;
    static const int coorddimension = cdim;

    typedef ct ctype;
    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    UnimplementedGeometry ( std::string message = "Geometry not implemented" )
      : message_( std::move( message ) )
    {}

    GeometryType type () const { DUNE_THROW( NotImplemented, message_ ); }
    bool affine () const { return false; }

    int corners () const { DUNE_THROW( NotImplemented, message_ ); }
    GlobalCoordinate corner ( int i ) const { DUNE_THROW( NotImplemented, message_ ); }
    GlobalCoordinate center () const { DUNE_THROW( NotImplemented, message_ ); }

    GlobalCoordinate global ( const LocalCoordinate &local ) const { DUNE_THROW( NotImplemented, message_ ); }
    LocalCoordinate local ( const GlobalCoordinate &global ) const { DUNE_THROW( NotImplemented, message_ ); }

    ctype integrationElement ( const LocalCoordinate &local ) const { DUNE_THROW( NotImplemented, message_ ); }
    ctype volume () const { DUNE_THROW( NotImplemented, message_ ); }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const { DUNE_THROW( NotImplemented, message_ ); }
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const { DUNE_THROW( NotImplemented, message_ ); }

  private:
    std::string message_;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_UNIMPLEMENTEDGEOMETRY_HH
