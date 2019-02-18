#!/bin/bash
# Helper script for the video NL-Bayes denoiser. Given
# a noisy video it computes the tvl1 optical flow and then
# calls the denoiser.

SEQ=$1 # noisy sequence path, in printf format, e.g. /my/sequence/frame-%02d.tif
FFR=$2 # first frame
LFR=$3 # last frame
SIG=$4 # noise standard dev.
OUT=$5 # output folder
PRM=$6 # denoiser parameters (between quotes, e.g. "-px1 5 -pt1 4 -verbose 0" 

mkdir -p $OUT

# we assume that the binaries are in the same folder as the script
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# error checking {{{1
for i in $(seq $FFR $LFR);
do
	file=$(printf $SEQ $i)
	if [ ! -f $file ]
	then
		echo ERROR: $file not found
		exit 1
	fi
done

# compute optical flow {{{1
$DIR/tvl1flow-seq.sh $SEQ $FFR $LFR $OUT

# run denoising {{{1
echo "$DIR/vnlbayes \
-i $SEQ -f $FFR -l $LFR -sigma $SIG \
-fof $OUT/tvl1_%03d_f.flo -bof $OUT/tvl1_%03d_b.flo \
-deno $OUT/deno_%03d.tif -bsic $OUT/bsic_%03d.tif $PRM"

$DIR/vnlbayes \
 -i $SEQ -f $FFR -l $LFR -sigma $SIG \
 -fof $OUT/tvl1_%03d_f.flo -bof $OUT/tvl1_%03d_b.flo \
 -deno $OUT/deno_%03d.tif -bsic $OUT/bsic_%03d.tif $PRM

# vim:set foldmethod=marker:
