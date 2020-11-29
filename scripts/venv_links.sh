#!/bin/bash
cd bin
for name in $(ls ../pdbtool/scripts/*.py)
do 
fname=$(basename $name)
fname=${fname%.*}
ln -s $name $fname
done
