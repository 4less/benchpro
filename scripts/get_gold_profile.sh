#!/bin/bash

GOLD=$1
OTU2ACC=$2
ACC2CLASS=$3
TMP=$4

GOLDTMP=$TMP/gold.tmp

mkdir $TMP

python scripts/join.py \
	-1 $GOLD \
	-2 $OTU2CLASS \
	-d "\t" \
	-j1 5 \
	-j2 0 \
	-c1 4,5 \
	-c2 1 > $GOLDTMP

python scripts/join.py \
	-1 $GOLDTMP \
	-2 $ACC2CLASS \
	-j1 2 \
	-j2 0 \
	-c1 0,1,2 \
	-c2 1 \
