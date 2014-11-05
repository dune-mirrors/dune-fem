dune_install(io.hh parameter.hh)

dune_add_subdirs(file streams visual matlab test)

dune_add_library(io OBJECT io.cc)
