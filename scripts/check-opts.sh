#!/bin/bash

# check for parameter pointing to DUNE base directory
# ---------------------------------------------------

DUNECONTROL="dune-common/bin/dunecontrol"

if test \( $# -lt 3 \) -o ! -e $1/$DUNECONTROL ; then
  echo "Usage: $0 <dune-base-dir> <opts-file> <headercheck{true,false}>"
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
OPTS=$2
HEADERCHECK=$3

errors=0

if ! test -e $OPTSDIR/$OPTS; then 
  echo "Error: $OPTS file not found in $OPTSDIR."
  exit 1
fi
#search for dune modules:

MODULES=""
for file in $(find -maxdepth 2 -name "dune.module" ) ; do
  MODULENAME="$(basename $(dirname $file))"
  MODULES+=" $MODULENAME"
done

# check for dependencies
MAKE_CHECK_DEPS=""
MAKE_CHECK_DEPS="$(source $OPTSDIR/$OPTS; echo $MAKE_CHECK_DEPS)"

MISSINGDEPS=""
for dep in $MAKE_CHECK_DEPS ; do
  if ! echo $MODULES | grep $dep >/dev/null; then
    MISSINGDEPS+=" $dep"
  fi
done

# obtain MAKE_FLAGS if specified
MAKE_CHECK_FLAGS=""
MAKE_CHECK_FLAGS="$(source $OPTSDIR/$OPTS; echo $MAKE_FLAGS)"

# we can skip this opts file if a dependend module is not found
if ! test -z "$MISSINGDEPS" ; then 
  echo "Skipping $OPTS due to missing dependencies:$MISSINGDEPS"
  exit 0
fi

BUILDDIR="$DUNEDIR/$(source $OPTSDIR/$OPTS; echo $BUILD_DIR)"

# perform tests
# -------------

do_configure()
{
  local check=`mktemp -p $WORKINGDIR check.XXXXXX`
  {
    $DUNECONTROL --opts=$OPTSDIR/$OPTS --builddir=$BUILDDIR all
    echo $? > $check
  } 2>&1 | dd conv=notrunc > $CONFIGLOG 2>/dev/null
  local return_value=`cat $check`
  rm $check
  return $return_value
}


echo
echo "Configuring with $OPTS ...."

CONFIGLOG="$WORKINGDIR/${OPTS%.opts}-conf.log"
if ! do_configure ; then
  echo "Error: Cannot configure with $OPTS (see ${OPTS%.opts}-conf.log)"
  exit 1
fi


# run headercheck
# ---------------

if $HEADERCHECK; then
  cd $WORKINGDIR
  if ! ./check-headers.sh $BUILDDIR dune-fem "$MAKE_CHECK_FLAGS"; then
    errors=$((errors+1))
  fi
fi
# run tests
# ---------

cd $WORKINGDIR
CHECKLOG="$WORKINGDIR/${OPTS%.opts}-check.log"
TESTLOG="$WORKINGDIR/${OPTS%.opts}-test.log"
TIMINGS="$WORKINGDIR/${OPTS%.opts}-timings.log"

if ! $SCRIPTSDIR/check-tests.sh $BUILDDIR/dune-fem "$MAKE_CHECK_FLAGS"; then
  build_test_errors=$(grep error: $WORKINGDIR/check-tests.log | wc -l)
  if test $build_test_errors -gt 0; then
    echo "Error: Building tests failed with $OPTS (see ${OPTS%.opts}-check.log)"
  fi

  if test -e $WORKINGDIR/LastTestsFailed.log; then
    failed_tests=$(cat $WORKINGDIR/LastTestsFailed.log | cut -d ':' -f2 | tr '\n' ' ' )
    echo "Error: For $OPTS the following tests failed: $failed_tests"
    echo "       see ${OPTS%.opts}-test.log for details."
    rm $WORKINGDIR/LastTestsFailed.log
  fi
  errors=$((errors+1))
fi
echo "Timings for $opts are stored in ${OPTS%.opts}-timings.log"
echo

# copy files
mv $WORKINGDIR/check-tests.log $CHECKLOG &>/dev/null
mv $WORKINGDIR/LastTest.log $TESTLOG &>/dev/null
mv $WORKINGDIR/CTestCostData.txt $TIMINGS &>/dev/null


if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
