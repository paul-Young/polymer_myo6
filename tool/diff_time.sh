#!/bin/bash

DATA_DIR=../output
DATA_FILE=diff_time.dat

cd $DATA_DIR
if [ -f $DATA_FILE ]
then 
	rm $DATA_FILE
fi

for file in $(ls)
do
	tail -n 2 $file | head -n 1 >> $DATA_FILE
done

mv $DATA_FILE ../Data
