SET(HEADERS
  attr.h
  cast.h
  common.h
  complex.h
  descr.h
  extensions.h
  eigen.h
  functional.h
  numpy.h
  operators.h
  pybind11.h
  pytypes.h
  stl.h
  stl_bind.h
  typeid.h
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/pypind11)
