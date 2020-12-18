dune_add_subdirs(test)

#install header files
dune_install( bindguard.hh
              coordinate.hh
              explicitfieldvector.hh
              fmatrixcol.hh
              forloop.hh
              geometrytypemap.hh
              hybrid.hh
              intersectionside.hh
              localcontribution.hh
              memory.hh
              referencevector.hh
              stackallocator.hh
              tupletypetraits.hh
              tupleutility.hh
              typeindexedtuple.hh
              utility.hh )

exclude_from_headercheck(
  referencevector.hh
)
