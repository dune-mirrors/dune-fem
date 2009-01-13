# if documentation was enabled 
if BUILD_DOCS
DOXYGEN_HTML = html 
DOXYGEN_DIST = html-dist

# what files does the doxygen-generated stuff depend on (find-syntax)
DEPENDON = -name \*.cc -o -name \*.hh -o -name \*.png
ADDITIONALDOCS = -name Doxyfile -o -name modules -o -name \*.dist -o -name \*.dev

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

# ----------------------------------------------------------------------------

doxygendir = $(datadir)/doc/dune-fem/doxygen

DISTCLEANFILES = doxygen.out doxygen-dist.out doxygen.log doxygen-dist.log
# we need this in the distribution
EXTRA_DIST = $(DOXYGEN_DIST) $(PAGES) Doxyfile doxygen-dist.out stamp
# these files should be removed on 'make clean'
CLEANFILES = doxygen.out doxygen-dist.out doxygen.log doxygen-dist.log

#all: $(DOXYGEN_HTML)
all: html
#	@echo "Please use 'make html' to build the documentation."

FORCE:

# check dependencies ourself to be robust (make stamp newer than any dependency)
stamp: FORCE
	@set -e; \
	if test ! -e stamp \
	|| test x"`find . \( $(ADDITIONALDOCS) \) -a -newer stamp`" != x \
	|| test x"`find $(top_srcdir) \( $(DEPENDON) \) -a -newer stamp`" != x ; \
	then \
	  touch stamp; \
	fi

html: doxygen.out

html-dist: doxygen-dist.out

doxygen.out: Doxyfile stamp
	@echo "Creating documentation (running doxygen)..."
	@( cat Doxyfile; $(DOXYFILE_OVERWRITE_HTML) ) > Doxyfile.tmp
	@$(DOXYGEN) Doxyfile.tmp > doxygen.log 2>&1
	@rm -f Doxyfile.tmp
	@cat doxygen.out

doxygen-dist.out: Doxyfile stamp
	@echo "Creating distribution documentation (running doxygen)..."
	@( cat Doxyfile; $(DOXYFILE_OVERWRITE_DIST) ) > Doxyfile.tmp
	@$(DOXYGEN) Doxyfile.tmp > doxygen-dist.log 2>&1
	@rm -f Doxyfile.tmp
	@cat doxygen-dist.out


# setting like in dune-web
BASEDIR=../..
CURDIR=doc/doxygen
EXTRAINSTALL = doxygen.log

# EXTRAINSTALL in the html subdirectory
DOXYGENINSTALL = *.html *.css *.png *.gif

install-data-local: $(DOXYGEN_DIST)
	$(mkinstalldirs) $(DESTDIR)/$(doxygendir); \
	list="$(srcdir)/html-dist/*.html $(srcdir)/html-dist/*.css $(srcdir)/html-dist/*.png $(srcdir)/html-dist/*.gif"; \
	for p in $$list; do \
		$(install_sh_DATA) $$p $(DESTDIR)/$(doxygendir); \
	done

uninstall-local:
	rm -f $(DESTDIR)/$(doxygendir)/*.html $(DESTDIR)/$(doxygendir)/*.css \
	$(DESTDIR)/$(doxygendir)/*.png $(DESTDIR)/$(doxygendir)/*.gif

# as all doxygen files fill be installed in the html subdirectory
# the BASEDIR has to be adjusted
%.build: FORCE
	$(MAKE) $* BASEDIR=../../..

# how to install the doxygen documentation
doxygen-install:
	if test -d $(DUNEWEBDIR) && test -r $(DUNEWEBDIR)/Make.global ; then \
	  pushd html; \
	    $(MAKE) -f $(DUNEWEBDIR)/Make.global \
          EXTRAINSTALL="$(DOXYGENINSTALL)" CURDIR="$(CURDIR)/dune-fem-html" install ; \
	  popd; \
	fi

web-install: doxygen-install

#include $(top_srcdir)/am/webstuff

clean-local:
	if test -e Doxydep; then rm -rf html html-dist stamp; fi
	rm -rf latex *~

# end BUILD_DOCS
endif

include $(top_srcdir)/am/global-rules
