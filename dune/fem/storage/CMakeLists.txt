dune_install(dynamicarray.hh
             eigenvector.hh
             envelope.hh
             objectstack.hh
             referencecounter.hh
             referencevector.hh
             rowreferencevector.hh
             singleton.hh
             singletonlist.hh
             subvector.hh)

dune_add_subdirs(test)

target_sources(dunefem PRIVATE singleton.cc)
