dune_install(alugridwriter.hh bartonnackmaninterface.hh boundaryidentifier.hh
             boundaryidprovider.hh capabilities.hh checkgeomaffinity.hh
             codimmap.hh commhelper.hh debug.hh
             double.hh engineconcept.hh femeoc.hh
             femeoctable.hh femtimer.hh fieldmatrixhelper.hh
             flops.hh fmatrixconverter.hh functor.hh
             griddeclaration.hh gridname.hh gridobjectstreams.hh
             gridsolution.hh gridwidth.hh h1norm.hh
             iteratorprovider.hh l1norm.hh l2norm.hh
             linesegmentsampler.hh lpnorm.hh mapgeomtype.hh
             metaprogramming.hh mpimanager.hh nil.hh
             nonconformitylevel.hh referenceelementprovider.hh subobjects.hh
             suite.hh test.hh threaditerator.hh
             threadmanager.hh validator.hh)

dune_add_subdirs(petsc threads)

dune_add_library(misc OBJECT femtimer.cc test.cc suite.cc)
