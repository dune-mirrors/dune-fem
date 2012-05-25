#!/bin/bash 
AUTHORS=$1
FILE=$2
OUT=tmp.cc.$LOGNAME 
echo "/**************************************************************************" > $OUT
echo " " >> $OUT 
cat $AUTHORS >> $OUT 
echo " " >> $OUT 
echo "**************************************************************************/" >> $OUT
cat $FILE >> $OUT
cp $OUT $FILE 
rm -f $OUT
