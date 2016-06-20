// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <string>
#include <tuple>

#include <dune/common/fmatrix.hh>
#include <dune/common/std/utility.hh>

#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/operators.h>

namespace py = pybind11;

template<class K, int rows, int cols>
static void registerFieldMatrix(py::handle scope)
{
  typedef Dune::FieldMatrix<K, rows, cols> FM;

  static const std::string clsName = "FieldMatrix" + std::to_string(rows) + std::to_string(cols);
  py::class_<FM> cls(scope, clsName.c_str());

  cls.def(py::init<>());

  cls.def("__init__",
      [] (FM& m, py::list l) {
        using std::size_t;

        new(&m) FM(K(0));

        size_t l_rows = l.size();
        for (size_t r = 0; r < l_rows; r++)
        {
          size_t l_cols = l[r].cast<py::list>().size();
          for (size_t c = 0; c < l_cols; c++)
            m[r][c] = l[r].cast<py::list>()[c].cast<K>();
        }
      });

  cls.def("__getitem__",
      [] (FM& m, std::size_t i) -> Dune::FieldVector<K, cols>& {
        if (i < rows)
          return m[i];
        else
          throw py::index_error();
      },
      py::return_value_policy::reference_internal);

  cls.def("__setitem__",
      [] (FM& m, std::size_t i, py::object l) {
        if(i < rows)
        {
          Dune::FieldVector<K, cols> v = l.cast<Dune::FieldVector<K, cols>>();
          std::size_t size = std::min(cols, (int) v.size()); // TODO: ask about relevance

          for (unsigned int j = 0; j < size; j++)
            m[i][j] = v[j];
        }
        else
          throw py::index_error();
      });

  cls.def("__len__", [] (const FM& m) -> std::size_t { return m.size(); });
  cls.def("invert", &FM::invert);

  cls.def(py::self += py::self);
  cls.def(py::self -= py::self);
  cls.def(py::self *= K());
  cls.def(py::self /= K());

  cls.def(py::self == py::self);
  cls.def(py::self != py::self);

  cls.def("__repr__",
      [] (const FM& m) {
        std::string repr = "DUNE FieldMatrix:\n(";

        for(int r = 0; r < rows; r++)
        {
          repr += "(";
          for (int c = 0; c < cols; c++)
            repr += (c > 0 ? ", " : "") + std::to_string(m[r][c]);
          repr += std::string(")") + (r < rows - 1 ? "\n" : "");
        }

        repr += ")";

        return repr;
      });

  cls.def_property_readonly("frobenius_norm",     [](const FM& m) { return m.frobenius_norm(); });
  cls.def_property_readonly("frobenius_norm2",    [](const FM& m) { return m.frobenius_norm2(); });
  cls.def_property_readonly("infinity_norm",      [](const FM& m) { return m.infinity_norm(); });
  cls.def_property_readonly("infinity_norm_real", [](const FM& m) { return m.infinity_norm_real(); });
  cls.def_property_readonly("rows", [](const FM& m) -> std::size_t { return rows; });
  cls.def_property_readonly("cols", [](const FM& m) -> std::size_t { return cols; });

  py::implicitly_convertible<py::list, FM>();
}

template<class K, int rows, int... cols>
static void registerCols(py::handle scope, std::integer_sequence<int, cols...>)
{
  std::ignore = std::make_tuple((registerFieldMatrix<K, rows, cols>(scope), 0)...);
}

template<class K, int... rows>
static void registerFieldMatrix(py::handle scope, std::integer_sequence<int, rows...> seq)
{
  std::ignore = std::make_tuple((registerCols<K, rows>(scope, seq), 0)...);
}
