#!/bin/bash

DATA_DIR=../Data
DATA_FILE=diff_time.dat

cd $DATA_DIR
if [ -f $DATA_FILE ]
then 
	rm $DATA_FILE
fi

for file in $(ls)
do
	tail -n 1 $file | awk '{print $1}' >> $DATA_FILE
done
