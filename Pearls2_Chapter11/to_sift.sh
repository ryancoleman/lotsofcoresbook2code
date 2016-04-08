#!/bin/bash
#
# ToSift.sh
# Created a script for extracting sift features from a set of images

BIN_PATH=/home/sanjana/visualsearch
IMG_DIR=/home/sanjana/visualsearch/pics

COUNTER=0
SIFT=$BIN_PATH/sift
OUTPUT_DIR=/home/sanjana/Desktop/IS/data

for f in $(find $IMG_DIR -name '*.pgm')
do
	$SIFT <$f >$OUTPUT_DIR/$COUNTER'.key'
	awk '$1~/^[0-9]+$/' $OUTPUT_DIR/$COUNTER'.key' | sed '1d' > 		$OUTPUT_DIR/$COUNTER
	COUNTER=$[$COUNTER +1]
done





