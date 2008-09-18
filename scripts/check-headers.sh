#!/bin/bash

PROJECTDIR=`pwd`
cd `dirname $0`
SCRIPTSDIR=`pwd`

if test ! -f $PROJECTDIR/dune.module ; then
  echo "This script can only be executed from a DUNE project directory."
  exit 1
fi

compile=1
for arg in $@ ; do
  if test "x$arg" == "xfast" ; then
    compile=0
    continue
  fi

  echo "Usage: $0 [fast]"
  exit 1
done

cd $PROJECTDIR
sourcedirs=`grep "^[[:space:]]*Sources:" dune.module | cut -s -d ':' -f2`
if test "x$sourcedirs" = "x" ; then
  echo -n "[.] "
  echo "Warning: 'Sources' not found in dune.module, assuming 'dune'."
  sourcedirs="dune"
fi

files=""
for directory in $sourcedirs ; do
  dirfiles=`find -H $directory -type d | sed "/$directory.*\/\..*/ d"`
  files="$files $dirfiles"
done

errors=0
for directory in $files ; do
  makefile=$directory/Makefile.am
  headerfiles=`ls -1 $directory | sed '/.*\.hh$/ p ; d'`

  if test ! -f $makefile ; then
    if test "x$headerfiles" != "x" ; then
      echo -n "[$directory] "
      echo "Warning: No 'Makefile.am' found."
    fi
    continue
  fi

  extradist=`cat $makefile | sed 'H ; s/.*//g ; x ; s/\n/ /g ; s/[ \t][ \t]*/ /g ; /\\\\$/! { p ; d } ; s/\\\\$// ; x ; d' | grep '^[ \t]*EXTRA_DIST' | sed 's/^.*EXTRA_DIST.*=//'`
  for pattern in $extradist ; do
    if test x`echo $pattern | grep '^[$]'` != x ; then
      echo -n "[$directory] "
      echo "Warning: Skipping check for variable '$pattern' in EXTRA_DIST."
      continue
    fi
    for file in $directory/$pattern ; do
      if test ! -f $file ; then
        echo -n "[$directory] "
        echo "Error: Distributed file does not exists: $file"
        errors=$((errors+1))
      fi
    done
  done

  ignore=`grep "IGNORE_HEADERS" $makefile`
  if test "x$ignore" != "x" ; then
    continue
  fi
 
  headers=`cat $makefile | sed 'H ; s/.*//g ; x ; s/\n/ /g ; s/[ \t][ \t]*/ /g ; /\\\\$/! { p ; d } ; s/\\\\$// ; x ; d' | grep '_HEADERS' | sed 's/^.*_HEADERS.*=//'`

  for header in $headers ; do
    if test x`echo $header | grep '^[$]'` != x ; then
      echo -n "[$directory] "
      echo "Warning: Skipping check for variable '$header' in headers declaration."
      continue
    fi
    if test ! -e $directory/$header ; then
      echo -n "[$directory] "
      echo "Error: Header does not exist: $header"
      errors=$((errors+1))
    fi
  done

  for filename in $headerfiles ; do
    found="false"
    for header in $headers ; do
      if test "x$header" = "x$filename" ; then
        found="true"
      fi
    done
    if test "x$found" = "xfalse" ; then
      echo -n "[$directory] "
      echo "Warning: Header is not included in Makefile.am: $filename"
    fi
  done

  if test $compile -ne 0 ; then
    for filename in $headerfiles ; do
      if ! $SHELL $SCRIPTSDIR/check-header.sh $directory/$filename &>/dev/null ; then
        echo -n "[$directory] "
        echo "Error: Header does not compile: $filename"
        errors=$((errors+1))
      fi
    done
  fi

#  for filename in `ls $directory | sed "/.*.hh/ p ; d"` ; do
#     echo $directory : $filename
#   done
done

if [ $errors -gt 0 ] ; then
  echo "Number of errors: $errors"
  exit 1
fi
