dune_install(io.hh parameter.hh)

dune_add_subdirs(file parameter streams visual test)

dune_add_library(io OBJECT io.cc)
