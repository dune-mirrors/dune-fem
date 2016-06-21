/*
    pybind11/extensions.h: Extensions to the C++11 python binding
    generator library for dune-fempy

    Copyright (c) 2016 Andreas Dedner <a.s.dedner@warwick.ac.uk>
    Copyright (c) 2016 Martin Nolte <nolte@mathematik.uni-freiburg.de>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "pybind11.h"

NAMESPACE_BEGIN(pybind11)

template <class T>
inline bool already_registered() {
  return static_cast<bool>(detail::get_type_info(typeid(T)));
}

NAMESPACE_END(pybind11)
