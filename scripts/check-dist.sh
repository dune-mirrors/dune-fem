#!/bin/bash

# check for parameter pointing to DUNE base directory
# ---------------------------------------------------

DUNECONTROL="dune-common/bin/dunecontrol"

if test \( $# -lt 1 \) -o ! -e $1/$DUNECONTROL ; then
  echo "Usage: $0 <dune-base-dir>"
  exit 1
fi

# set up some variables
# ---------------------

WORKINGDIR=`pwd`
cd $1
DUNEDIR=`pwd`
FEMDIR="$DUNEDIR/dune-fem"
SCRIPTSDIR="$FEMDIR/scripts"
OPTSDIR="$SCRIPTSDIR/opts"

MODULES="dune-common dune-grid dune-istl dune-fem"

# fetch missing tarballs from website
# -----------------------------------

for MODULE in $MODULES ; do
  cd $DUNEDIR/$MODULE
  if test x`find -maxdepth 1 -name "*.tar.gz"` != x ; then
    continue
  fi

  if test "$MODULE" != "dune-fem" ; then
    echo
    echo "Downloading tarball for $MODULE from dune-project.org..."
    wget -q "http://www.dune-project.org/download/1.0beta7/$MODULE-1.0beta7.tar.gz"
  fi

  if test x`find -maxdepth 1 -name "*.tar.gz"` == x ; then
    echo "Fatal: No tarball available for $MODULE"
    exit 1
  fi
done

# perform tests
# -------------

TESTDIR=`mktemp -d -p $WORKINGDIR dune-tmp-XXXXXX`

errors=0
for OPTS in `cd $OPTSDIR ; ls *.opts` ; do
  echo "Checking $OPTS..."

  cd $TESTDIR
  rm -rf *
 
  find $DUNEDIR -maxdepth 2 -name "*.tar.gz" -exec tar -xzf \{\} \;

  for MODULE in $MODULES ; do
    mv `ls -d ${MODULE}*` $MODULE
  done

  CONFIGLOG="$WORKINGDIR/${OPTS%.opts}-conf.out"
  if ! $DUNECONTROL --opts=$OPTSDIR/$OPTS all &> $CONFIGLOG ; then
    echo "Error: Cannot configure with $OPTS (see $CONFIGLOG)"
    errors=$((errors+1))
    continue
  fi

  CONFIGOUT="$WORKINGDIR/${OPTS%.opts}-check.out"
  cd dune-fem/fem/test
  if ! make check &> $CHECKLOG ; then
    echo "Error: Check failed with $OPTS (see $CHECKLOG)"
    errors=$((errors+1))
  fi
done

# clean up
# --------

rm -rf $TESTDIR

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
