# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = test

streamsdir =  $(includedir)/dune/fem/io/streams
streams_HEADERS = asciistreams.hh \
                  streams.hh \
                  streams_inline.hh \
		  virtualstreams.hh \
		  xdrstreams.hh

include $(top_srcdir)/am/global-rules

