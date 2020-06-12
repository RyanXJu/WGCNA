#!/bin/bash

# module load R/3.6.1

DIR=$(dirname "$0")

# help info
Help(){
  echo ""
  echo "Usage: ./WGCNA_part1.sh -o [Output Path] -d [Expression data] -t [Traits data] -l [log base] -c [add constant before log]"
  echo ""
	echo "--------------------- $(tput bold)WGCNA_part1$(tput sgr0) -------------------------------"
	echo "This part will: "
	echo "1. remove genes with too many missing values, or low variance"
	echo "2. build sample tree with kept genes for user to cut outliers in the following part."
	echo ""
	echo "--------------------- $(tput bold)Data$(tput sgr0)---------------------------------------"
	echo "$(tput bold)Input data$(tput sgr0)"
	echo "Expression data      normalized gene expression data (csv/tsv file, rows--genes, columns--samples)"
	echo "Trait data           trait data of all the samples (csv/tsv file, rows--samples, columns--trait)"
	echo ""
  echo ""
	echo "$(tput bold)Output data$(tput sgr0)"
	echo "sampleTree.Rdata                        Data needed for WGCNA_part2.sh"
	echo "Part1_SampleClustering.pdf              Visualization of sample tree"
	echo ""
	echo "--------------------- $(tput bold)Parameters$(tput sgr0) --------------------------------" 
  echo "-o|--output       directory of output"
	echo "-d|--data         path of expression data"
	echo "-t|--trait        path of trait data"
	echo "-l|--log          base of Log transformation, no log transformation if not defined"
	echo "-c|--constant     add constant to expression values before log tranformation"
  echo""
}


# parameters
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
        Help
        exit 1
        ;;
    *)    # unknown option
        shift # past argument
        ;;
    esac
done


# check for needed parameters
if [ -z "$DATA" ] ; then
	echo -e "$(tput setaf 1)NO expression data supplied!\n$(tput sgr 0)Please enter the path of expression data (i.e. -d path/datExpr.tsv)"
	exit 1
fi

if [ -z "$TRAIT" ] ; then
	echo -e "$(tput setaf 1)NO trait data supplied!\n$(tput sgr 0)Please enter the path of trait data (i.e. -t path/datTraits.tsv)"
	exit 1
fi

if [ -z "$OUTPUT" ] ; then
	echo -e "$(tput setaf 1)NO output path supplied!\n$(tput sgr 0)Please define output directory (i.e. -o output_path)" 
	exit 1
fi

# run R
Rscript "${DIR}/WGCNA_part1.R" $OUTPUT $DATA $TRAIT $LOG $CONSTANT
