#!/bin/bash

if [ $# -ne 2 ]; then
    echo second
    echo $0: usage: run-peak-detector.sh input_folder output_folder
    exit 1
fi

in_folder=$1
out_folder=$2

# creates the out folder
rm -r $out_folder
mkdir $out_folder

#copy files
cp $in_folder/* $out_folder

# runs the stuff
for f in $out_folder/*
do
  ./peax $f $f
done

