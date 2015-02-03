dune_install(adaptivefunction.hh blockvectorfunction.hh
             combinedfunction.hh productfunction.hh vectorfunction.hh
             subfunction.hh blockvectordiscretefunction.hh petscdiscretefunction.hh)

dune_add_subdirs(common localfunction adaptivefunction
                 combinedfunction blockvectorfunction
                 productfunction vectorfunction blockvectordiscretefunction
                 blockvectors petscdiscretefunction)
