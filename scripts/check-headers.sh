#!/bin/bash

WORKINGDIR=`pwd`
cd `dirname $0`
SCRIPTSDIR=`pwd`
cd "$SCRIPTSDIR/.."
FEMDIR=`pwd`

compile=1
for arg in $@ ; do
  if test "x$arg" == "xfast" ; then
    compile=0
    continue
  fi

  echo "Usage: $0 [fast]"
  exit 1
done

errors=0
for directory in `find fem -type d | sed "/fem.*\/\..*/ d"` ; do
  makefile=$directory/Makefile.am

  if test -f $makefile ; then
    extradist=`cat $makefile | sed 'H ; s/.*//g ; x ; s/\n/ /g ; s/[ \t][ \t]*/ /g ; /\\\\$/! { p ; d } ; s/\\\\$// ; x ; d' | grep '^.*EXTRA_DIST' | sed 's/^.*EXTRA_DIST.*=//'`
    for pattern in $extradist ; do
      for file in $directory/$pattern ; do
        if test ! -f $file ; then
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
      if test ! -e $directory/$header ; then
        echo "Error: Header does not exist: $directory/$header"
        errors=$((errors+1))
      fi
    done

    for filename in `ls $directory | sed '/.*.hh/ p ; d'` ; do
      found="false"
      for header in $headers ; do
        if test "x$header" = "x$filename" ; then
          found="true"
        fi
      done
      if test "x$found" = "xfalse" ; then
        echo "Warning: Header is not included in Makefile.am: $directory/$filename"
      fi
    done

    if test $compile -ne 0 ; then
      for filename in `ls $directory | sed '/.*.hh/ p ; d'` ; do
        cd $WORKINGDIR
        if ! $SHELL $SCRIPTSDIR/check-header.sh $directory/$filename &>/dev/null ; then
          echo "Error: Header does not compile: $directory/$filename"
          errors=$((errors+1))
        fi
      done
      cd $FEMDIR
    fi

  else
    echo "Warning: Makefile.am not found in $directory"
  fi

#  for filename in `ls $directory | sed "/.*.hh/ p ; d"` ; do
#     echo $directory : $filename
#   done
done

if [ $errors -gt 0 ] ; then
  echo "Number of errors: $errors"
  exit 1
fi
