dune_install(alugridwriter.hh bartonnackmaninterface.hh
             boundaryidprovider.hh capabilities.hh checkgeomaffinity.hh
             compatibility.hh debug.hh double.hh femeoc.hh
             femeoctable.hh femtimer.hh fieldmatrixhelper.hh
             flops.hh fmatrixconverter.hh functor.hh
             griddeclaration.hh gridname.hh gridobjectstreams.hh
             gridsolution.hh gridwidth.hh h1norm.hh
             l1norm.hh l2norm.hh gridfunctionview.hh
             linesegmentsampler.hh lpnorm.hh mapgeomtype.hh
             metaprogramming.hh mpimanager.hh
             nonconformitylevel.hh umfpack.hh domainintegral.hh)

dune_add_subdirs(petsc threads)

target_sources(dunefem PRIVATE femtimer.cc)

exclude_from_headercheck( compatibility.hh )
