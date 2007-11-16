#!/bin/bash

# check for parameter pointing to DUNE base directory
# ---------------------------------------------------

DUNECONTROL="dune-common/bin/dunecontrol"

if test \( $# -lt 1 \) -o ! -e $1/$DUNECONTROL ; then
  echo "Usage: $0 <dune-base-dir>"
  exit 1
fi

WORKINGDIR=`pwd`
cd $1
DUNEDIR=`pwd`

# check headers in Makefile.am
# ----------------------------

# configure with minimal options
# ------------------------------

OPTSDIR="$DUNEDIR/dune-fem/scripts/opts"
MINIMALOPTS="$OPTSDIR/minimal.opts"

if test ! -e $MINIMALOPTS ; then
  echo "Error: $MINIMALOPTS not found."
  exit 1
fi

cd $DUNEDIR
if ! $DUNECONTROL --opts=$MINIMALOPTS all &> $WORKINGDIR/minimal-svn-conf.out ; then
  echo "Error: Cannot configure with minimal options (see $WORKINGDIR/minimal-svn-conf.out)."
  exit 1
fi

# build tarballs
# --------------

MODULES="dune-common dune-grid dune-istl dune-fem"

for MODULE in $MODULES ; do
  echo "Making tarball in $MODULE..."

  cd $DUNEDIR/$MODULE
  find -maxdepth 1 -name "*.tar.gz" -delete
  if ! make dist &> $WORKINGDIR/$MODULE-dist.out ; then
    echo "Error: Cannot make tarball for $MODULE (see $WORKINGDIR/$MODULE-dist.out)"
    exit 1
  fi
done

# check documentation
# -------------------

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
fi
echo "Done."
exit 0

#--------------------------------------------------------------------

FILENAME=$1
DUNEDIR=$2
FEMDIR="$DUNEDIR/dune-fem"
SCRITPSDIR="$FEMDIR/scripts"
MYDIR=$3

cd $SCRITPSDIR
# file to mail
if test -e $FILENAME ; then
  rm -f $FILENAME
fi
touch $FILENAME

cd $SCRITPSDIR
if ! ./test-preconf.sh $FEMDIR ; then 
  echo "*******************" >> $FILENAME
  echo "PRECONFIGURE TESTS:" >> $FILENAME
  echo >> $FILENAME  
  cat preconf.out >> $FILENAME
  echo "*******************" >> $FILENAME  
  echo >> $FILENAME  
fi

# hier noch einfuegen: for Schleife ueber config files
cd $SCRITPSDIR
CONFIGFILE="$SCRIPTSDIR/opts/config_all.opts"
if ! ./test-conf.sh $CONFIGFILE $FEMDIR ; then 
  echo "*******************" >> $FILENAME
  echo "CONFIGURE TESTS:" >> $FILENAME
  echo >> $FILENAME  
  cat conf.out >> $FILENAME
  echo "*******************" >> $FILENAME  
  echo >> $FILENAME  
fi

cd $SCRITPSDIR
if ! ./test-postconf.sh $FEMDIR ; then 
  echo "*******************" >> $FILENAME
  echo "POSTCONFIGURE TESTS:" >> $FILENAME
  echo >> $FILENAME  
  cat postconf.out >> $FILENAME
  echo "*******************" >> $FILENAME  
  echo >> $FILENAME  
fi

cd $SCRITPSDIR
if ! ./test-run.sh $FEMDIR ; then 
  echo "*******************" >> $FILENAME
  echo "RUN TESTS:" >> $FILENAME
  echo >> $FILENAME  
  cat run.out >> $FILENAME
  echo "*******************" >> $FILENAME  
  echo >> $FILENAME  
fi

cd $SCRITPSDIR
mv $FILENAME $MYDIR/
exit 0


