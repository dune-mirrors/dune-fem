// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_FEMPY_GEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_FEMPY_GEOMETRY_REFERENCEELEMENTS_HH

// #include <experimental/array>
// 4.9 doesn't have experimental array - reimplemented here:
namespace std
{
  template<class T, class... Tail, class Elem = typename std::decay<T>::type>
  std::array<Elem,1+sizeof...(Tail)> make_array(T&& head, Tail&&... values)
  {
    return { std::forward<T>(head), std::forward<Tail>(values)... };
  }
}


#include <dune/geometry/referenceelements.hh>

#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{
  namespace FemPy
  {

    namespace py = pybind11;

    // registerReferenceElement
    // ------------------------

    template<class ctype, int dim>
    auto registerReferenceElement_(py::handle scope)
    {
      typedef typename Dune::ReferenceElement<ctype, dim> RefElement;

      static const std::string name = "ReferenceElement" + std::to_string(dim);
      py::class_<RefElement> cls(scope, name.c_str());

      cls.def("size", [](const RefElement& e, int c)                { return e.size(c); });
      cls.def("size", [](const RefElement& e, int i, int c, int cc) { return e.size(i, c, cc); });
      cls.def("type", [](const RefElement& e, int i, int c)         { return e.type(i, c); });

      cls.def("subEntity"             , &RefElement::subEntity);
      cls.def("position"              , &RefElement::position);
      cls.def_property_readonly("center", [](const RefElement &ref) { return ref.position(0,0); } );
      // cls.def("checkInside"           , &RefElement::checkInside);
      cls.def("volume"                , &RefElement::volume);
      cls.def("integrationOuterNormal", &RefElement::integrationOuterNormal);

      cls.def_property_readonly("type", [](const RefElement& e) { return e.type(); });

      // registerGeometry<RefElement>(cls, std::make_integer_sequence<int, dim+1>());

      return cls;
    }

    template<class ctype, int... dim>
    void registerReferenceElement(py::handle scope, std::integer_sequence<int, dim...>)
    {
      std::ignore = std::make_tuple(registerReferenceElement_<ctype, dim>(scope)...);
    }



    // registerReferenceElements
    // -------------------------

    template<class ctype, int dim>
    auto registerReferenceElements_(py::handle scope)
    {
      typedef typename Dune::ReferenceElements<ctype, dim> RefElements;

      static const std::string name = "ReferenceElements" + std::to_string(dim);
      py::class_<RefElements> cls(scope, name.c_str());

      cls.def_static("general", &RefElements::general, py::return_value_policy::reference_internal);
      cls.def_static("simplex", &RefElements::simplex, py::return_value_policy::reference_internal);
      cls.def_static("cube"   , &RefElements::cube,    py::return_value_policy::reference_internal);

      return py::cast(RefElements());
    }

    template<class ctype, int... dim>
    void registerReferenceElements(py::module module, std::integer_sequence<int, dim...>)
    {
      static const auto referenceElements =
        std::make_array(registerReferenceElements_<ctype, dim>(module)...);

      module.def("ReferenceElements", [](int dimension) { return referenceElements[dimension]; });
    }

  } // namespace FemPy

} // namespace Dune

#endif // ifndef DUNE_FEMPY_REFERENCEELEMENTS_HH
