# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

storagedir =  $(includedir)/dune/fem/storage
storage_HEADERS = dynamicarray.hh \
                  eigenvector.hh \
                  envelope.hh \
                  objectstack.hh \
                  referencecounter.hh \
                  singletonlist.hh \
                  subvector.hh \
                  eigenvector.hh

SUBDIRS = test

include $(top_srcdir)/am/global-rules
