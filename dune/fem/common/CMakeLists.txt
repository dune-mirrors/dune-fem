#install header files
dune_install( coordinate.hh
              fmatrixcol.hh
              referencevector.hh
              stackallocator.hh
              tupleforeach.hh
              tupletypetraits.hh
              tupleutility.hh
              typeindexedtuple.hh
              utility.hh )

dune_add_subdirs(test)
