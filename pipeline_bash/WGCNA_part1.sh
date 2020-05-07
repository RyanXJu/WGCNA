#!/bin/bash

module load R/3.6.1

LOG=""
CONSTANT=""

while [[ $# -gt 0 ]]
do
    key="${1}"
    case ${key} in
    -o|--output)
	OUTPUT="${2}"
	shift
	shift
	;;
    -d|--data)
        DATA="${2}"
        shift # past argument
        shift # past value
        ;;
    -t|--trait)
        TRAIT="${2}"
        shift # past argument
        shift # past value
        ;;
    -l|--log)
	LOG="${2}"
	shift
	shift
	;;
    -c|--constant)
	CONSTANT="${2}"
	shift
	shift
	;;
    -h|--help)
        echo "-o|--output : directory of output"
	echo "-d|--data : path of expression data"
	echo "-t|--trait : path of trait data"
	echo "-l|--log : base of Log transformation, no transformation will be performed if not define"
	echo "-c|--constant : constant add to value before log tranformation"
        shift # past argument
        ;;
    *)    # unknown option
        shift # past argument
        ;;
    esac
done

#echo $OUTPUT
#echo $DATA
#echo $TRAIT
#echo $LOG
#echo $CONSTANT


if [ -z "$DATA" ] ; then
	echo "NO expression data supplied!"
	exit 1
fi

if [ -z "$TRAIT" ] ; then
	echo "NO trait data supplied!"
	exit 1
fi

if [ -z "$OUTPUT" ] ; then
	echo "Please define output directory" 
	exit 1
fi


Rscript WGCNA_part1.R $OUTPUT $DATA $TRAIT $LOG $CONSTANT
