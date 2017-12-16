#!/bin/bash
# Downloads and extracts the data from the given URL 1
in_folder="input"
out_folder="output"

rm -r $in_folder
mkdir $in_folder

rm candy_*.zip
wget http://www.imada.sdu.dk/~jbaumbac/download/teaching/ws17-18/DM847/project/toy_data/candy_raw.zip
unzip candy_*.zip -d $in_folder


# creates the out folder
rm -r $out_folder
mkdir $out_folder

#copy files
cp $in_folder/* $out_folder

# runs the stuff
for f in $out_folder/BD*.csv
do
  ./peax $f $f
done

