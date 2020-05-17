#ifndef DUNE_FEMPY_PYBIND11_PYBIND11_HH
#define DUNE_FEMPY_PYBIND11_PYBIND11_HH

#include <dune/fempy/pybind11/gridfunction.hh>
#include <dune/fempy/pybind11/space.hh>

#include <pybind11/complex.h>
#if HAVE_EIGEN
#include <pybind11/eigen.h>
#endif
#include <dune/python/extensions.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#endif // #ifndef DUNE_FEMPY_PYBIND11_PYBIND11_HH
