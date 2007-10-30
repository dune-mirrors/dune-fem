# if documentation was enabled 
#if BUILD_DOCS
DOXYGENTAG = doxygen-tag 
WHAT = # bestpractice.html
DOXYGEN_HTML = html 
DOXYGEN_DIST = html-dist

DISTCLEANFILES = doxygen.log doxygen-tag doxygen.out

all: $(WHAT) $(DOXYGEN_HTML)

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/doxygen
EXTRAINSTALL = doxygen.log

# EXTRAINSTALL in the html subdirectory
DOXYGENINSTALL = *.html *.css *.png *.gif

if WML
# themeing doxygen
DOXYGENHEADER = doxy-header.html
DOXYGENHEADERBUILD = $(DOXYGENHEADER).build
DOXYGENFOOTER = doxy-footer.html
DOXYGENSTYLESHEET = dune-doxy.css
DOYXFILE_OVERWRITE = echo "HTML_HEADER = $(DOXYGENHEADER)"; \
  echo "HTML_FOOTER = $(DOXYGENFOOTER)"; \
  echo "HTML_STYLESHEET = $(DOXYGENSTYLESHEET)";
DOXYFILESHORT_OVERWRITE = $(DOYXFILE_OVERWRITE) \
  echo "SHORT_NAMES = YES";
endif


#we need this in the distribution
EXTRA_DIST = $(WHAT) $(DOXYGEN_DIST) $(DOXYGENHEADER) doxy-footer.html $(PAGES)

doxygendir = $(datadir)/doc/dune-fem/doxygen
doxygen_DATA = $(WHAT)


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
ADDITIONALDOCS = -name Doxyfile -o -name docrules -o -name mainpage -o -name modules -o -name svnpage -o -name suggestions 

# as all doxygen files fill be installed in the html subdirectory
# the BASEDIR has to be adjusted
%.build: FORCE
	$(MAKE) $* BASEDIR=../../..

# check dependency ourself to be robust
doxygen-tag: FORCE $(DOXYGENHEADERBUILD)
	set -e; \
	if test -e Doxyfile; then \
	  if test ! -e doxygen-tag || test x"`find $(top_srcdir) \( $(DEPENDON) \) -a -newer doxygen-tag -print`" != x || test x"`find $(BASEDIR) \( $(ADDITIONALDOCS) \) -a -newer doxygen-tag -print`" != x ; then \
	    echo Running doxygen. This may take a while... ; \
	    ( cat Doxyfile; $(DOYXFILE_OVERWRITE) ) > Doxyfile.tmp; \
	    $(DOXYGEN) Doxyfile.tmp > doxygen.log 2>&1 <&- ; \
	    rm -f Doxyfile.tmp; \
	    touch doxygen-tag ; \
	  fi; \
	fi

html: doxygen-tag

html-dist: $(DOXYGENHEADERBUILD)
	rm -rf html
	set -e; \
	if test -e Doxyfile; then \
	  ( cat Doxyfile; $(DOXYFILESHORT_OVERWRITE) ) > Doxyfile.tmp; \
	  echo Running doxygen. This may take a while... ; \
	  $(DOXYGEN) Doxyfile.tmp > doxygen.log 2>&1 <&- ; \
	  rm -f Doxyfile.tmp; \
	fi
	if test -d html; then rm -rf html-dist; mv html html-dist; fi

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

include $(top_srcdir)/am/webstuff

clean-local:
	if test -e Doxydep; then rm -rf html html-dist; fi
	rm -rf latex *~

CLEANFILES = $(WHAT) $(DOXYGENTAG) $(DOXYGENHEADER)

# end BUILD_DOCS
#endif

include $(top_srcdir)/am/global-rules
