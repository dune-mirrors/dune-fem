#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <dune-fem-dir>"
  exit 1
fi

cd $1

errors=0
for directory in `find fem -type d | sed "/fem.*\/\..*/ d"` ; do
  makefile=$directory/Makefile.am

  if test -f $makefile ; then
    ignore=`grep "IGNORE_HEADERS" $makefile`

    if test "x$ignore" != "x" ; then
      #echo "Continuing: $directory"
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
