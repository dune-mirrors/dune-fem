#!/bin/bash 

KEY1=$1 
KEY2=$2
FILE=$3

echo "replacing '$KEY1' with '$KEY2' in $FILE"
CMD="'s|$KEY1|$KEY2|'"
#echo "sed $CMD"
# copy file first
cp $FILE $FILE.tmp
# pipe command to sh 
echo "cat $FILE.tmp | sed $CMD > $FILE" | sh 
rm $FILE.tmp
