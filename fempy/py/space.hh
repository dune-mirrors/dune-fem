#ifndef DUNE_FEMPY_PY_SPACE_HH
#define DUNE_FEMPY_PY_SPACE_HH

#include <dune/fem/misc/l2norm.hh>

#include <dune/corepy/common/fmatrix.hh>
#include <dune/corepy/common/fvector.hh>

#include <dune/corepy/pybind11/complex.h>
#include <dune/corepy/pybind11/pybind11.h>

#include <dune/fempy/function/virtualizedgridfunction.hh>

namespace Dune
{

  namespace FemPy
  {

    template <class GridPart, class Field, int dimR>
    double distance
    (Dune::FemPy::VirtualizedGridFunction<GridPart, Dune::FieldVector<Field,dimR>> &u,
     Dune::FemPy::VirtualizedGridFunction<GridPart, Dune::FieldVector<Field,dimR>> &v)
    {
      Dune::Fem::L2Norm<GridPart> norm( u.gridPart(), 8);
      return norm.distance(u,v);
    }

    // registerSpace
    // -------------

    namespace detail
    {
      template< class Space, class Cls >
      void registerSpace ( pybind11::module module, Cls &cls )
      {
        typedef typename Space::GridPartType GridPartType;
        typedef typename Space::FunctionSpaceType::RangeFieldType RangeFieldType;
        static const int dimRange = Space::dimRange;
        if (!std::is_same<RangeFieldType,double>::value)
        {
          registerFieldVector<RangeFieldType>(module, std::make_integer_sequence<int, 10>());
          registerFieldMatrix<RangeFieldType>(module, std::make_integer_sequence<int, 5>());
        }

        typedef typename Space::GridPartType GridPart;

        cls.def_property_readonly( "grid", [](Space &sp) -> const GridPart& {return sp.gridPart();} );

        cls.def( "__init__", [] ( Space &instance, GridPart &grid ) {
            new( &instance ) Space( grid );
          }, pybind11::keep_alive< 1, 2 >() );
        cls.def( "distance", [] ( Space&,
          Dune::FemPy::VirtualizedGridFunction<GridPartType, Dune::FieldVector<RangeFieldType,dimRange>> &u,
          Dune::FemPy::VirtualizedGridFunction<GridPartType, Dune::FieldVector<RangeFieldType,dimRange>> &v)
          {
            return distance(u,v);
          } );
      }
    }
    template< class Space, class Holder, class AliasType >
    void registerSpace ( pybind11::module module, pybind11::class_<Space,Holder,AliasType> &cls )
    {
      detail::registerSpace<Space>(module,cls);
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
