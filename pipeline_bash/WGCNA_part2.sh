#!/bin/bash

module load R/3.6.1

CUTHEIGHT=100000000
while [[ $# -gt 0 ]]
do
    key="${1}"
    case ${key} in
    -o|--output)
	OUTPUT="${2}"
	shift
	shift
	;;
	    -c|--cut)
        CUTHEIGHT="${2}"
        shift # past argument
        shift # past value
        ;;
    -h|--help)
        echo "-o|--output : directory of output"
        shift # past argument
        ;;
    *)    # unknown option
        shift # past argument
        ;;
    esac
done

if [ -z "$OUTPUT" ] ; then
	echo "Need the directory of sampleTree.RData!"
	exit 1
fi

echo $CUTHEIGHT

Rscript WGCNA_part2.R $OUTPUT $CUTHEIGHT
