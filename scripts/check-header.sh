#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <header>"
  exit 1
fi

# use new version of header check 
make HEADER=$1 headercheck 1>/dev/null
SUCCESS=$?

exit $SUCCESS 

#############################################
# old version of header check
#############################################

PROJECTDIR=`pwd`
cd `dirname $0`
SCRIPTSDIR=`pwd`
HEADER=$1

if test ! -f $PROJECTDIR/dune.module ; then
  echo "This script can only be executed from a DUNE project directory."
  exit 1
fi

if ! test -e $PROJECTDIR/$HEADER ; then
  echo "'$HEADER' does not exist."
  exit 1
fi

if test ! -f $PROJECTDIR/make.inc ; then
  echo "'make.inc' missing from project directory. Did you forget to configure?"
  exit 1
fi

TEMPDIR=`mktemp -d --tmpdir tmp-header-XXXXXX`

cd $TEMPDIR
cp $PROJECTDIR/make.inc Makefile
echo ".cc.o:" >> Makefile
echo -e -n "\t\${CXX} -c -I$PROJECTDIR " >> Makefile
echo '${AM_CPPFLAGS} ${AM_CXXFLAGS} -o $@ $<' >> Makefile

echo "#include <config.h>" >> test.cc
echo "#include <${HEADER}>" >> test.cc
echo "#include <${HEADER}>" >> test.cc
echo "int main () {}" >> test.cc

make test.o 1>/dev/null
SUCCESS=$?

rm -f test.o
rm -f test.cc
rm -f Makefile

cd $PROJECTDIR
rmdir $TEMPDIR

exit $SUCCESS
