dune_install(
  adaptcallbackhandle.hh
  adaptmanager.hh
  adaptationmanager.hh
  allgeomtypes.hh
  auxiliarydofs.hh
  basesetlocalkeystorage.hh
  cachedcommmanager.hh
  capabilities.hh
  commindexmap.hh
  commoperations.hh
  communicationmanager.hh
  datacollector.hh
  dataprojection.hh
  defaultcommhandler.hh
  discretefunctionspace.hh
  dofmanager.hh
  dofstorage.hh
  functionspace.hh
  functionspaceinterface.hh
  interpolate.hh
  loadbalancer.hh
  localinterpolation.hh
  localrestrictprolong.hh
  restrictprolongfunction.hh
  restrictprolonginterface.hh
  restrictprolongtuple.hh
  slavedofs.hh
  uniquefacetorientation.hh
)

dune_add_subdirs( dataprojection )

# avoid deprecation warning in header-check.
exclude_from_headercheck(
  adaptmanager.hh
  )
