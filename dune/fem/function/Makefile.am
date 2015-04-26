# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = common \
          localfunction \
          adaptivefunction \
          combinedfunction \
          blockvectorfunction \
          vectorfunction \
				  blockvectordiscretefunction \
				  blockvectors \
					petscdiscretefunction \
          test

functiondir =  $(includedir)/dune/fem/function
function_HEADERS = adaptivefunction.hh \
                   blockvectorfunction.hh \
                   combinedfunction.hh \
                   vectorfunction.hh \
                   subfunction.hh \
								   blockvectordiscretefunction.hh \
									 petscdiscretefunction.hh 

include $(top_srcdir)/am/global-rules

