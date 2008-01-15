# if documentation was enabled 
if BUILD_DOCS
DOXYGEN_HTML = html 
DOXYGEN_DIST = html-dist

DISTCLEANFILES = doxygen.log doxygen-tag doxygen-dist-tag doxygen.out doxygen-dist.out

all: $(DOXYGEN_HTML)

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/doxygen
EXTRAINSTALL = doxygen.log

# EXTRAINSTALL in the html subdirectory
DOXYGENINSTALL = *.html *.css *.png *.gif


DOXYFILE_OVERWRITE = echo "PROJECT_NAME = $(PACKAGE_NAME)"; \
                     echo "PROJECT_NUMBER = $(PACKAGE_VERSION)";
	
DOXYFILE_OVERWRITE_HTML = $(DOXYFILE_OVERWRITE) \
                          echo "HTML_OUTPUT = $(DOXYGEN_HTML)"; \
			  echo "WARN_LOGFILE = doxygen.out"; \
			  echo "FILE_PATTERNS += *.dev";
DOXYFILE_OVERWRITE_DIST = $(DOXYFILE_OVERWRITE) \
                          echo "HTML_OUTPUT = $(DOXYGEN_DIST)"; \
			  echo "WARN_LOGFILE = doxygen-dist.out"; \
			  echo "INTERNAL_DOCS = NO";

#we need this in the distribution
EXTRA_DIST = $(DOXYGEN_DIST) $(PAGES)

doxygendir = $(datadir)/doc/dune-fem/doxygen


install-data-local: $(DOXYGEN_DIST)
	$(mkinstalldirs) $(DESTDIR)/$(doxygendir); \
	list="$(srcdir)/html-dist/*.html $(srcdir)/html-dist/*.css $(srcdir)/html-dist/*.png $(srcdir)/html-dist/*.gif"; \
	for p in $$list; do \
		$(install_sh_DATA) $$p $(DESTDIR)/$(doxygendir); \
	done

uninstall-local:
	rm -f $(DESTDIR)/$(doxygendir)/*.html $(DESTDIR)/$(doxygendir)/*.css \
	$(DESTDIR)/$(doxygendir)/*.png $(DESTDIR)/$(doxygendir)/*.gif

# what files does the doxygen-generated stuff depend on (find-syntax)
DEPENDON = -name \*.cc -o -name \*.hh -o -name \*.png
ADDITIONALDOCS = -name Doxyfile -o -name modules -o -name \*.dist

# as all doxygen files fill be installed in the html subdirectory
# the BASEDIR has to be adjusted
%.build: FORCE
	$(MAKE) $* BASEDIR=../../..

# check dependency ourself to be robust
doxygen-tag: FORCE $(DOXYGENHEADERBUILD)
	@set -e; \
	if test -e Doxyfile; then \
	  if test ! -e doxygen-tag || test x"`find $(top_srcdir) \( $(DEPENDON) \) -a -newer doxygen-tag -print`" != x || test x"`find . \( $(ADDITIONALDOCS) -o -name \*.dev \) -a -newer doxygen-tag -print`" != x ; then \
	    echo "Creating documentation (running doxygen) ..." ; \
	    ( cat Doxyfile; $(DOXYFILE_OVERWRITE_HTML) ) > Doxyfile.tmp; \
	    $(DOXYGEN) Doxyfile.tmp > doxygen.log 2>&1 <&- ; \
	    rm -f Doxyfile.tmp; \
	    touch doxygen-tag ; \
	    cat doxygen.out ; \
	  fi; \
	fi


doxygen-dist-tag: FORCE $(DOXYGENHEADERBUILD)
	@set -e; \
	if test -e Doxyfile; then \
	  if test ! -e doxygen-dist-tag || test x"`find $(top_srcdir) \( $(DEPENDON) \) -a -newer doxygen-dist-tag -print`" != x || test x"`find . \( $(ADDITIONALDOCS) \) -a -newer doxygen-dist-tag -print`" != x ; then \
	    echo "Creating distribution documentation (running doxygen) ..." ; \
	    ( cat Doxyfile; $(DOXYFILE_OVERWRITE_DIST) ) > Doxyfile.tmp; \
	    $(DOXYGEN) Doxyfile.tmp > doxygen.log 2>&1 <&- ; \
	    rm -f Doxyfile.tmp; \
	    touch doxygen-dist-tag ; \
	    cat doxygen-dist.out ; \
	  fi \
	fi

html: doxygen-tag

html-dist: doxygen-dist-tag

# how to install the doxygen documentation
doxygen-install:
	if test -d $(DUNEWEBDIR) && test -r $(DUNEWEBDIR)/Make.global ; then \
	  pushd html; \
	    $(MAKE) -f $(DUNEWEBDIR)/Make.global \
          EXTRAINSTALL="$(DOXYGENINSTALL)" CURDIR="$(CURDIR)/dune-fem-html" install ; \
	  popd; \
	fi

web-install: doxygen-install

# nice trick from the GNU make infopage to force a rule to run
FORCE:

#include $(top_srcdir)/am/webstuff

clean-local:
	if test -e Doxydep; then rm -rf html html-dist; fi
	rm -rf latex *~

CLEANFILES = doxygen-tag doxygen-dist-tag

# end BUILD_DOCS
endif

include $(top_srcdir)/am/global-rules
