#ifndef DUNE_FEMPY_PY_GRID_GEOMETRY_HH
#define DUNE_FEMPY_PY_GRID_GEOMETRY_HH

#include <string>

#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // PyCornerRange
    // -------------

    template< class Geometry >
    struct PyCorners
    {
      PyCorners ( const Geometry &geometry, pybind11::object ref )
        : geometry_( geometry ), ref_( ref )
      {}

      const Geometry &geometry () { return geometry_; }

    private:
      const Geometry &geometry_;
      pybind11::object ref_;
    };


    // PyCornerIterator
    // ----------------

    template< class Geometry >
    struct PyCornerIterator
    {
      PyCornerIterator ( const PyCorners< Geometry > corners ) : corners_( corners ) {}

      typename Geometry::GlobalCoordinate next ()
      {
        if( index_ == corners_.geometry().corners() )
          throw pybind11::stop_iteration();
        return corners_.geometry().corner( index_++ );
      }

    private:
      PyCorners< Geometry > corners_;
      int index_ = 0;
    };



    // registerGridGeometry
    // --------------------

    template< class Geometry >
    pybind11::class_< Geometry > registerGridGeometry ( pybind11::handle scope )
    {
      static const std::string clsName = "Geometry" + std::to_string( Geometry::mydimension );
      pybind11::class_< Geometry > cls( scope, clsName.c_str() );

      pybind11::class_< PyCornerIterator< Geometry > > itCls( cls, "CornerIterator" );
      itCls.def( "__iter__", [] ( PyCornerIterator< Geometry > &it ) -> PyCornerIterator< Geometry > & { return it; } );
      itCls.def( "__next__", &PyCornerIterator< Geometry >::next );

      pybind11::class_< PyCorners< Geometry > > cCls( cls, "Corners" );
      cCls.def( "__iter__", [] ( const PyCorners< Geometry > &c ) { return PyCornerIterator< Geometry >( c ); } );

      cls.def_property_readonly( "corners", [] ( pybind11::object geo ) {
          return PyCorners< Geometry >( geo.cast< const Geometry & >(), geo );
        } );

      cls.def_property_readonly( "center", &Geometry::center );
      cls.def_property_readonly( "volume", &Geometry::volume );

      cls.def_property_readonly( "affine", &Geometry::affine );

      cls.def( "position", &Geometry::global );
      cls.def( "integrationElement", &Geometry::integrationElement );

      cls.def_property_readonly( "domain", [] ( const Geometry &geo ) -> const Dune::ReferenceElement< typename Geometry::ctype, Geometry::mydimension > & {
          return Dune::ReferenceElements< typename Geometry::ctype, Geometry::mydimension >::general( geo.type() );
        }, pybind11::return_value_policy::reference );

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GEOMETRY_HH
