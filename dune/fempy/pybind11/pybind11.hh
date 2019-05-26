#ifndef DUNE_FEMPY_PYBIND11_PYBIND11_HH
#define DUNE_FEMPY_PYBIND11_PYBIND11_HH

#include <dune/fempy/pybind11/gridfunction.hh>
#include <dune/fempy/pybind11/space.hh>

#include <dune/python/pybind11/complex.h>
#if HAVE_EIGEN
#include <dune/python/pybind11/eigen.h>
#endif
#include <dune/python/pybind11/extensions.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#endif // #ifndef DUNE_FEMPY_PYBIND11_PYBIND11_HH
