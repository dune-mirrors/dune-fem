# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = common \
          localfunction \
	  adaptivefunction \
	  attachedfunction \
	  combinedfunction \
	  blockvectorfunction \
	  productfunction \
	  vectorfunction

functiondir =  $(includedir)/dune/fem/function
function_HEADERS = adaptivefunction.hh \
                   attachedfunction.hh \
                   blockvectorfunction.hh \
		   combinedfunction.hh \
		   productfunction.hh \
		   vectorfunction.hh

include $(top_srcdir)/am/global-rules

