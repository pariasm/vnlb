#!/bin/bash
# Computes forward and backward optical flow for a sequence

SEQ=$1          # sequence path
FFR=$2          # first frame
LFR=$3          # last frame
OUT=$4          # output folder
DIR=${5:-"all"} # direction of flow

mkdir -p $OUT

# we assume that the binaries are in the same folder as the script
BINDIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

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
TVL1="$BINDIR/tvl1flow"
FSCALE=1

if [ "$DIR" == "bwd" ] || [ "$DIR" == "all" ]
then
	for i in $(seq $((FFR+1)) $LFR);
	do
		file=$(printf $OUT"/tvl1_%03d_b.flo" $i)
		if [ ! -f $file ]
		then
			$TVL1 $(printf $SEQ $i) \
					$(printf $SEQ $((i-1))) \
					$file \
					0 0.25 0.2 0.3 100 $FSCALE 0.5 5 0.01 0; 
		fi
	done
	cp $(printf $OUT"/tvl1_%03d_b.flo" $((FFR+1))) $(printf $OUT"/tvl1_%03d_b.flo" $FFR)
fi

if [ "$DIR" == "fwd" ] || [ "$DIR" == "all" ]
then
	for i in $(seq $FFR $((LFR-1)));
	do
		file=$(printf $OUT"/tvl1_%03d_f.flo" $i)
		if [ ! -f $file ]
		then
			$TVL1 $(printf $SEQ $i) \
					$(printf $SEQ $((i+1))) \
					$file \
					0 0.25 0.2 0.3 100 $FSCALE 0.5 5 0.01 0; 
		fi
	done
	cp $(printf $OUT"/tvl1_%03d_f.flo" $((LFR-1))) $(printf $OUT"/tvl1_%03d_f.flo" $LFR)
fi

