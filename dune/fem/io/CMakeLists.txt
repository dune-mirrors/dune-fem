dune_install(io.hh parameter.hh)

dune_add_subdirs(file parameter streams visual test)

target_sources(dunefem PRIVATE io.cc)
