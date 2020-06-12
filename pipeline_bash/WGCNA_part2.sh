#!/bin/bash

# module load R/3.6.1

DIR=$(dirname "$0")

Help(){
  echo ""
  echo "Usage: ./WGCNA_part2.sh -o [PATH of sampleTree.Rdata] -c [Height]"
  echo ""
	echo "--------------------- $(tput bold)WGCNA_part2$(tput sgr0) -------------------------------"
	echo "This part will: "
	echo "1. Cut sample tree to remove outliers according to user's decision"
	echo "2. Prepare scale-free topology information for user to choose soft-threshold in the following part."
	echo ""
	echo "--------------------- $(tput bold)Data$(tput sgr0)---------------------------------------"
	echo "$(tput bold)Input data$(tput sgr0)"
	echo "sampleTree.Rdata (returned by WGCNA_part1)"
	echo ""
	echo "$(tput bold)Output data$(tput sgr0)"
	echo "topology.Rdata                          Data needed for WGCNA_part3.sh"
	echo "Part2_SampleClustering_cut.pdf          Visualization of sample tree"
	echo "Part2_SampleDendrogram.pdf              Visualization of the kept samples and their traits"
	echo "Part2_ScaleFreeTopology.pdf             Topplogy and connectivity of the network"
	echo ""
	echo "--------------------- $(tput bold)Parameters$(tput sgr0) --------------------------------" 
	echo "-o|--output      directory of output (folder of sampleTree.Rdata)"
	echo "-c|--cutHeight   height to cut sample tree. The largest cluster will be kept for analysis. (use all samples if not define) "
  echo ""
}

CUTHEIGHT=""
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
        Help
        exit 1
        ;;
    *)    # unknown option
        shift # past argument
        ;;
    esac
done

# check parameters
if [ -z "$OUTPUT" ] ; then
	echo -e "$(tput setaf 1)NO output path supplied!\n$(tput sgr 0)Please enter the directory of sampleTree.RData! (i.e. -o path)"
	exit 1
fi

if [ ! -z "$CUTHEIGHT" ] && [[ ! "$CUTHEIGHT" =~ ^[0-9]+$ ]]  ; then  # CUTHEIGHT is number 
	echo -e "$(tput setaf 1)Cut height needs to be a positive number!\n$(tput sgr 0)Please consult Part1_SampleClustering.pdf!"
	exit 1
fi


Rscript "${DIR}/WGCNA_part2.R" $OUTPUT $CUTHEIGHT
