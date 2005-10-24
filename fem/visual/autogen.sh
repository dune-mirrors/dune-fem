#!/bin/sh

set -e

# may be used to force a certain automake-version e.g. 1.7
AMVERS=

if test x$1 = "x" ; then
  echo "Usage: ./autogen.sh DUNEDIR"
  exit 1
fi

if test ! -d $1/m4 ; then
  echo $1/m4 not found! Wrong directory supplied?
  exit 1
fi

if test "x$AMVERS" != x ; then
  echo Warning: explicitly using automake version $AMVERS
  # binaries are called automake-$AMVERS
  AMVERS="-$AMVERS"
fi

# convert to absolute path so that aclocal 1.8 does the right thing
DUNEM4=`cd $1/m4 && pwd`

aclocal$AMVERS -I $DUNEM4 -I ./

libtoolize --force

autoheader

automake$AMVERS --add-missing

autoconf
