miscdir =  $(includedir)/dune/fem/misc

SUBDIRS = petsc test

misc_HEADERS = alugridwriter.hh \
               bartonnackmaninterface.hh \
               boundaryidentifier.hh \
               boundaryidprovider.hh \
               capabilities.hh \
               checkgeomaffinity.hh \
               codimmap.hh \
               combineinterface.hh \
               commhelper.hh \
               debug.hh \
               double.hh \
               engineconcept.hh \
               femeoc.hh \
               femeoctable.hh \
               field.hh \
               functor.hh \
               femtimer.hh \
               femtuples.hh \
               fieldmatrixhelper.hh \
               griddeclaration.hh \
               gridsolution.hh \
               gridwidth.hh \
               h1norm.hh \
               l1error.hh \
               l2error.hh \
               l1norm.hh \
               l2norm.hh \
               lpnorm.hh \
               linesegmentsampler.hh \
               mapgeomtype.hh \
               metaprogramming.hh \
               mpimanager.hh \
               nonconformitylevel.hh \
               referenceelementprovider.hh \
               selectvalue.hh \
               subobjects.hh \
               suite.cc \
               suite.hh \
               test.cc \
               test.hh \
               utility.hh \
               validator.hh \
               fmatrixconverter.hh \
               gridname.hh \
               threaditerator.hh \
               threadmanager.hh

noinst_LTLIBRARIES = libmisc.la
libmisc_la_SOURCES = femtimer.cc test.cc suite.cc

include $(top_srcdir)/am/global-rules
