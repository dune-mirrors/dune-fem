dune_add_subdirs(local)

dune_install(hdivprojection.hh dgl2projection.hh
             vtxprojection.hh)

# avoid deprecation warnings during headercheck
exclude_from_headercheck(
  hdivprojection.hh
  dgl2projection.hh
  l2projection.hh
  )
