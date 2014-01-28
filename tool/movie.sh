#!/bin/bash

# first input should run number, write a check routine later

script_dir=$(pwd)
struct_dir=../Struct_data/$1
movie_name=sim$1.pdb

cd $struct_dir

# in struct directory

echo "MODEL" > $movie_name
for struct in $(ls|grep -v $movie_name)
do

	cat $struct >> $movie_name
	echo "ENDMDL" >> $movie_name
	echo "MODEL">> $movie_name
	 
done
echo "ENDMDL" >> $movie_name

mv $movie_name $script_dir
	

