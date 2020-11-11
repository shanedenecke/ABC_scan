#!/bin/bash
# Determine the current encoding of the file
encoding=$(file -i "$2" | sed "s/.*charset=\(.*\)$/\1/")
 
if [ ! "$1" == "${encoding}" ]
then
# Encodings differ, we have to encode
echo "recoding from ${encoding} to $1 file : $2"
recode ${encoding} $1 $2
fi
