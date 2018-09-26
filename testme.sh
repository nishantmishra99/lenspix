#!/bin/bash
count=$1

file=input.txt
if [ -f $file ] ; then
 rm $file
fi
echo $count
for i in $(seq $count)
do
echo test.py $i >> $file
done

