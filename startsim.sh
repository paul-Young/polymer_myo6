#!/bin/bash

BIN=run_myo

if [ ! -f $BIN ] 
then
	make
fi

if [ ! -d output ] 
then
	mkdir output 
fi


for run in {1..4}
do
	nice ./run_myo pdbmyoVI1.pdb pdbmyoVI1.pdb $run > output/out$run.log &
done

