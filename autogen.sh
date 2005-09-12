#!/bin/sh
# $Id$

#### barf on errors
set -e

# may be used to force a certain automake-version e.g. 1.7
AMVERS=

# everybody who checks out the CVS wants the maintainer-mode to be enabled
# (should be off for source distributions, this should happen automatically)
#
DEFAULTCONFOPT="--enable-maintainer-mode"

# default values
DEBUG=1
OPTIM=0

usage () {
    echo "Usage: ./autogen.sh [options]"
    echo "  -i, --intel        use intel compiler"
    echo "  -g, --gnu          use gnu compiler (default)"
    echo "  --opts=FILE        use compiler-options from FILE"
    echo "  -d, --debug        switch debug-opts on"
    echo "  -n, --nodebug      switch debug-opts off"
    echo "  -o, --optim        switch optimization on"
    echo "  --with-dune=PATH   directory with dune/ inside"
    echo "  -h, --help         you already found this :)"
    echo
    echo "Parameters not in the list above are directly passed to configure. See"
    echo
    echo "    ./configure --help"
    echo
    echo "for a list of additional options"
}

# no compiler set yet
COMPSET=0
for OPT in $* ; do

    set +e
    # stolen from configure...
    # when no option is set, this returns an error code
    arg=`expr "x$OPT" : 'x[^=]*=\(.*\)'`
    set -e

    case "$OPT" in
	-i|--intel)   . ./icc.opts ; COMPSET=1 ;;
	-g|--gnu)     . ./gcc.opts ; COMPSET=1 ;;
	--mpi)     . ./mpi.opts    ; COMPSET=1 ;;
	--opts=*)
	    if [ -r $arg ] ; then
	      echo "reading options from $arg..."
	      . ./$arg ;
	      COMPSET=1;
	    else
	      echo "Cannot open compiler options file $arg!" ;
	      exit 1;
	    fi ;;
	-d|--debug)   DEBUG=1 ;;
	-n|--nodebug) DEBUG=0 ;;
	-o|--optim)   OPTIM=1 ;;
	-h|--help) usage ; exit 0 ;;
	# special hack: use the with-dune-dir for aclocal-includes
	--with-dune=*)
	    eval DUNEDIR=$arg
	    # add the option anyway
	    CONFOPT="$CONFOPT $OPT" ;;
	# pass unknown opts to ./configure
	*) CONFOPT="$CONFOPT $OPT" ;;
    esac
done

# set special m4-path if --with-dune is set
if [ x$DUNEDIR != x ] ; then
    # aclocal from automake 1.8 seems to need an absolute path for inclusion
    FULLDIR=`cd $DUNEDIR && pwd`

    # automagically use directory above if complete Dune-dir was supplied
    if test `basename $FULLDIR` = "dune" ; then
      FULLDIR=`cd $FULLDIR/.. && pwd`
    fi

    ACLOCALOPT="-I $FULLDIR/dune/m4/ -I ./" 
fi

# use the free compiler as default :-)
if [ "$COMPSET" != "1" ] ; then
    echo "No compiler set, using GNU compiler as default"
    . ./gcc.opts
fi

# create flags
COMPFLAGS="$FLAGS"

# maybe add debug flag
if [ "$DEBUG" = "1" ] ; then	
    COMPFLAGS="$COMPFLAGS $DEBUGFLAGS"
fi

# maybe add optimization flag
if [ "$OPTIM" = "1" ] ; then	
    COMPFLAGS="$COMPFLAGS $OPTIMFLAGS"
fi

# check if automake-version was set
if test "x$AMVERS" != x ; then
  echo Warning: explicitly using automake version $AMVERS
  # binaries are called automake-$AMVERS
  AMVERS="-$AMVERS"
fi

#### create all autotools-files

echo "--> libtoolize..."
# force to write new versions of files, otherwise upgrading libtools
# doesn't do anything...
libtoolize --force

echo "--> aclocal..."
aclocal$AMVERS $ACLOCALOPT

# sanity check to catch missing --with-dune
if ! grep DUNE aclocal.m4 > /dev/null ; then
    echo "aclocal.m4 doesn't contain any DUNE-macros, this would crash autoconf"
    echo "or automake later. Maybe you should provide a --with-dune=PATH parameter"
    exit 1
fi

echo "--> autoheader..."
autoheader

echo "--> automake..."
automake$AMVERS --add-missing

echo "--> autoconf..."
autoconf

#### start configure with special environment

export CC="$COMP"
export CXX="$CXXCOMP"
export CPP="$COMP -E"

export CFLAGS="$COMPFLAGS"
export CXXFLAGS="$COMPFLAGS"

./configure $DEFAULTCONFOPT $CONFOPT

