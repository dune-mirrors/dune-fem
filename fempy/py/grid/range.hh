#error DON'T INCLUDE THIS - TAKEN CATE OF BY COREPY

#ifndef DUNE_FEMPY_PY_GRID_RANGE_HH
#define DUNE_FEMPY_PY_GRID_RANGE_HH

#include <string>
#include <utility>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // PyGridPartRange
    // ---------------

    template< class GridPart, int codim >
    struct PyGridPartRange
    {
      typedef typename GridPart::template Codim< codim >::IteratorType Iterator;

      PyGridPartRange ( const GridPart &gridPart, pybind11::object ref )
        : gridPart_( gridPart ), ref_( std::move( ref ) )
      {}

      Iterator begin () const { return gridPart_.template begin< codim >(); }
      Iterator end () const { return gridPart_.template end< codim >(); }

    private:
      const GridPart &gridPart_;
      pybind11::object ref_;
    };



    // PyGridPartIterator
    // ------------------

    template< class GridPart, int codim >
    struct PyGridPartIterator
    {
      typedef PyGridPartRange< GridPart, codim > Range;
      typedef typename GridPart::template Codim< codim >::EntityType Entity;

      PyGridPartIterator ( const Range &range ) : range_( range ), it_( range_.begin() ) {}

      Entity next ()
      {
        if( it_ == range_.end() )
          throw pybind11::stop_iteration();

        Entity entity = *it_;
        ++it_;
        return entity;
      }

    private:
      Range range_;
      typename Range::Iterator it_;
    };



    // registerPyGridPartRange
    // -----------------------

    template< class GridPart, int codim >
    void registerPyGridPartRange ( pybind11::handle scope, const char *rgName )
    {
      typedef PyGridPartRange< GridPart, codim > Range;
      typedef PyGridPartIterator< GridPart, codim > Iterator;

      static const std::string itName = std::string( rgName ) + "Iterator";
      pybind11::class_< Iterator > itCls( scope, itName.c_str() );
      itCls.def( "__iter__", [] ( Iterator &it ) -> Iterator & { return it; } );
      itCls.def( "__next__", &Iterator::next );

      pybind11::class_< Range > rgCls( scope, rgName );
      rgCls.def( "__iter__", [] ( const Range &range ) { return Iterator( range ); } );
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_RANGE_HH
