dune_install(asciistreams.hh binarystreams.hh sionlibstreams.hh
             standardstreams.hh streams.hh streams_inline.hh
             tuples.hh virtualstreams.hh xdrstreams.hh)

dune_add_subdirs(test)

exclude_from_headercheck(
  xdrstreams.hh
  )
