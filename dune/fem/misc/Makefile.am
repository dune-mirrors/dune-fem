miscdir =  $(includedir)/dune/fem/misc

SUBDIRS = petsc test threads

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
	       femtimer.cc \
	       femtimer.hh \
	       field.hh \
	       fieldmatrixhelper.hh \
	       fmatrixconverter.hh \
	       functor.hh \
	       griddeclaration.hh \
	       gridname.hh \
	       gridobjectstreams.hh \
	       gridsolution.hh \
	       gridwidth.hh \
	       h1norm.hh \
				 iteratorprovider.hh \
	       l1error.hh \
	       l1norm.hh \
	       l2error.hh \
	       l2norm.hh \
	       linesegmentsampler.hh \
	       lpnorm.hh \
	       mapgeomtype.hh \
	       metaprogramming.hh \
	       mpimanager.hh \
	       nil.hh \
	       nonconformitylevel.hh \
	       referenceelementprovider.hh \
	       subobjects.hh \
	       suite.cc \
	       suite.hh \
	       test.cc \
	       test.hh \
	       threaditerator.hh \
	       threadmanager.hh \
	       validator.hh

noinst_LTLIBRARIES = libmisc.la
libmisc_la_SOURCES = femtimer.cc test.cc suite.cc

include $(top_srcdir)/am/global-rules
