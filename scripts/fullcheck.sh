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

echo
echo "Host Name: $HOSTNAME"
echo "Host Type: $HOSTTYPE"

# set up some variables
# ---------------------

WORKINGDIR=`pwd`
cd $1
DUNEDIR=`pwd`
FEMDIR="$DUNEDIR/dune-fem"
SCRIPTSDIR="$FEMDIR/scripts"
OPTSDIR="$SCRIPTSDIR/opts"

MODULES="dune-common dune-grid dune-istl dune-fem"

errors=0

# check headers in Makefile.am
# ----------------------------

echo
echo "Checking Makefile.am's *_HEADERS variables..."
cd $FEMDIR
if ! $SCRIPTSDIR/check-headers.sh fast ; then
  errors=$((errors+1))
fi

# configure with minimal options
# ------------------------------

MINIMALOPTS="$OPTSDIR/minimal.opts"

if test ! -e $MINIMALOPTS ; then
  echo "Error: $MINIMALOPTS not found."
  exit 1
fi

minimal_configure()
{
  local check=`mktemp -p $WORKINGDIR check.XXXXXX`
  {
    $DUNECONTROL --opts=$MINIMALOPTS all
    echo $? > $check
  } 2>&1 | dd conv=notrunc > $WORKINGDIR/minimal-svn-conf.out 2>/dev/null
  local return_value=`cat $check`
  rm $check
  return $return_value
}

echo
echo "Configuring with minimal options..."
cd $DUNEDIR
if ! minimal_configure ; then
  echo "Fatal: Cannot configure with minimal options (see $WORKINGDIR/minimal-svn-conf.out)."
  exit 1
fi


# check headers
# -------------

echo
echo "Checking headers..."
cd $FEMDIR
if ! $SCRIPTSDIR/check-headers.sh ; then
  errors=$((errors+1))
fi

# check documentation
# -------------------

echo
echo "Checking documentation..."
cd $WORKINGDIR
if ! $SCRIPTSDIR/check-doxygen.sh $FEMDIR ; then
  errors=$((errors+1))
fi

# build tarballs
# --------------

for MODULE in $MODULES ; do
  echo
  echo "Making tarball in $MODULE..."

  cd $DUNEDIR/$MODULE
  find -maxdepth 1 -name "*.tar.gz" -delete
  if ! make dist &> $WORKINGDIR/$MODULE-dist.out ; then
    echo "Error: Cannot make tarball for $MODULE (see $WORKINGDIR/$MODULE-dist.out)"
    if test $MODULE == dune-fem ; then
      errors=$((errors+1))
    fi
  fi
done

# check distributions
# -------------------

cd $WORKINGDIR
if ! $SCRIPTSDIR/check-dist.sh $DUNEDIR ; then
  errors=$((errors+1))
fi

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
