#!/bin/bash

# check for parameter pointing to DUNE base directory
# ---------------------------------------------------

DUNECONTROL="dune-common/bin/dunecontrol"

if test \( $# -lt 1 \) -o ! -e $1/$DUNECONTROL ; then
  echo "Usage: $0 <dune-base-dir>"
  exit 1
fi

echo "Full Check of dune-fem"
echo "----------------------"

# set up some variables
# ---------------------

WORKINGDIR=`pwd`
cd $1
DUNEDIR=`pwd`
FEMDIR="$DUNEDIR/dune-fem"
SCRIPTSDIR="$FEMDIR/scripts"
OPTSDIR="$SCRIPTSDIR/opts"

errors=0

# check headers in Makefile.am
# ----------------------------

echo
echo "Checking Makefile.am's *_HEADERS variables..."
if ! $SCRIPTSDIR/check-headers.sh $FEMDIR ; then
  errors=$((errors+1))
fi

# configure with minimal options
# ------------------------------

MINIMALOPTS="$OPTSDIR/minimal.opts"

if test ! -e $MINIMALOPTS ; then
  echo "Error: $MINIMALOPTS not found."
  exit 1
fi

echo
echo "Configuring with minimal options..."
cd $DUNEDIR
if ! $DUNECONTROL --opts=$MINIMALOPTS all &> $WORKINGDIR/minimal-svn-conf.out ; then
  echo "Error: Cannot configure with minimal options (see $WORKINGDIR/minimal-svn-conf.out)."
  exit 1
fi

# check documentation
# -------------------

echo
echo "Checking documentation..."
if ! $SCRIPTSDIR/check-doxygen.sh $FEMDIR ; then
  errors=$((errors+1))
fi

# build tarballs
# --------------

MODULES="dune-common dune-grid dune-istl dune-fem"

for MODULE in $MODULES ; do
  echo
  echo "Making tarball in $MODULE..."

  cd $DUNEDIR/$MODULE
  find -maxdepth 1 -name "*.tar.gz" -delete
  if ! make dist &> $WORKINGDIR/$MODULE-dist.out ; then
    errors=$((errors+1))
    echo "Error: Cannot make tarball for $MODULE (see $WORKINGDIR/$MODULE-dist.out)"
    exit 1
  fi
done

exit 0

# perform test builds
# -------------------

TESTDIR=`mktemp -p $WORKINGDIR dune-tmp-XXXXXX`
mkdir $TESTDIR

errors=0
for OPTS in `cd $OPTSDIR ; ls *.opts` ; do
  echo "Checking $OPTS..."

  cd $TESTDIR
  rm -rf *
 
  find $DUNEDIR -maxdepth 2 -name "*.tar.gz" -exec tar -xzf \{\} \;

  CONFIGOUT="$WORKINGDIR/${OPTS%.opts}-conf.out"
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

rm -rf $TESTDIR

# show results
# ------------

if test $errors -ne 0 ; then
  echo "Done ($errors occurred)."
  exit 1
else
  echo "Done."
  exit 0
fi
