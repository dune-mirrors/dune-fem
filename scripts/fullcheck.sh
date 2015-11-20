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

errors=0

NONOPTSFILE=""
if test -e $OPTSDIR/config_none.opts; then
  NONOPTSFILE="config_none.opts"
fi

OPTFILES=$(cd $OPTSDIR ; ls *.opts | grep -v $NONOPTSFILE)

for OPTS in $NONOPTSFILE $OPTFILES ; do
  echo
  echo "Performing checks for $OPTS..."

  cd $WORKINGDIR
  case $OPTS in
    $NONOPTSFILE)
      # perform headercheck only ones
      if ! $SCRIPTSDIR/check-opts.sh $DUNEDIR $OPTS true; then
        errors=$((errors+1))
      fi
    ;;
    *)
      # perform headercheck only ones
      if ! $SCRIPTSDIR/check-opts.sh $DUNEDIR $OPTS false; then
        errors=$((errors+1))
      fi
    ;;
    esac
done

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
