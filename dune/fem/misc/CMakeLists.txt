dune_install(alugridwriter.hh bartonnackmaninterface.hh
             boundaryidprovider.hh capabilities.hh checkgeomaffinity.hh
             codimmap.hh commhelper.hh debug.hh
             double.hh femeoc.hh
             femeoctable.hh femtimer.hh fieldmatrixhelper.hh
             flops.hh fmatrixconverter.hh functor.hh
             griddeclaration.hh gridname.hh gridobjectstreams.hh
             gridsolution.hh gridwidth.hh h1norm.hh
             iteratorprovider.hh l1norm.hh l2norm.hh
             linesegmentsampler.hh lpnorm.hh mapgeomtype.hh
             metaprogramming.hh mpimanager.hh
             nonconformitylevel.hh subobjects.hh
             validator.hh)

dune_add_subdirs(petsc threads)

dune_add_library(misc OBJECT femtimer.cc)
