SUBDIRS = pages

EXTRA_DIST = Doxyfile modules
CLEANFILES = stamp doxygen.out doxygen.log

# ----------------------------------------------------------------------------

# what files does the doxygen-generated stuff depend on (find-syntax)
DEPENDON = -name \*.cc -o -name \*.hh -o -name \*.png
ADDITIONALDOCS = -name Doxyfile -o -name modules -o -name pages\*

# Things to add to Doxyfile
DOXYFILE_OVERWRITE = echo "PROJECT_NAME = $(PACKAGE_NAME)"; \
                     echo "PROJECT_NUMBER = $(PACKAGE_VERSION)";

doxygendir = $(datadir)/doc/dune-fem/doxygen

# ----------------------------------------------------------------------------

if BUILD_DOCS
all: html/index.html
else
all:
endif

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

html/index.html: Doxyfile stamp
	@echo "Creating documentation (running doxygen)..."
	@( cat Doxyfile; $(DOXYFILE_OVERWRITE) ) > Doxyfile.tmp
	@$(DOXYGEN) Doxyfile.tmp > doxygen.log 2>&1
	@rm -f Doxyfile.tmp
	@cat doxygen.out

if BUILD_DOCS
install-data-local: html/index.html
	@echo "Installing documentation..."
	@set -e; \
	$(mkinstalldirs) $(DESTDIR)/$(doxygendir); \
	list="$(srcdir)/html/*.html $(srcdir)/html/*.css $(srcdir)/html/*.png $(srcdir)/html/*.gif"; \
	for p in $$list; do \
	  $(install_sh_DATA) $$p $(DESTDIR)/$(doxygendir); \
	done

uninstall-local:
	@echo "Uninstalling documentation..."
	@rm -f $(DESTDIR)/$(doxygendir)/*.html
	@rm -f $(DESTDIR)/$(doxygendir)/*.css
	@rm -f $(DESTDIR)/$(doxygendir)/*.png
	@rm -f $(DESTDIR)/$(doxygendir)/*.gif
endif

clean-local:
	rm -rf html

include $(top_srcdir)/am/global-rules
