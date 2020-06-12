#!/bin/bash

# module load R/3.6.1

DIR=$(dirname "$0")

Help(){
  echo ""
  echo "Usage: ./WGCNA_part3.sh -o [PATH of topology.Rdata] -s [Soft-threshold] -n [Network type]"
  echo ""
	echo "--------------------- $(tput bold)WGCNA_part3$(tput sgr0) -------------------------------"
	echo "This part will: "
	echo "1. Build scale-free network with chosen soft-threshold"
	echo "2. Detect co-expression moodules"
	echo "2. Calculate module-traits correlation"
	echo ""
	echo "--------------------- $(tput bold)Data$(tput sgr0)---------------------------------------"
	echo "$(tput bold)Input data$(tput sgr0)"
	echo "topology.Rdata (returned by WGCNA_part2)"
	echo ""
	echo "$(tput bold)Output data$(tput sgr0)"
	echo "network.RData                           Data needed for WGCNA_part4.sh"
	echo "WGCNA_expressionData.tsv                Expression data of all genes in the network"
	echo "WGCNA_traitData.tsv                     Trait data of all the kept samples"
	echo "WGCNA_moduleEigengenes.tsv              Module eigengenes expression in kept samples"
	echo "Part3_ScaleFreeTopology_sft.pdf         Network topplogy and connectivity level with the selected soft-threshold"
	echo "Part3_GeneClusterDendrogram.pdf         Visualization of gene clusters"
	echo "Part3_ModuleTraitCorrelation.pdf        Correlation heatmap of module eigengenes and traits"

	echo ""
	echo "--------------------- $(tput bold)Parameters$(tput sgr0) --------------------------------" 
	echo "-o|--output      directory of output (folder of topology.Rdata)"
	echo "-s|--sft         soft-threshold to build scale-free network "
	echo "-n|--network     OPTION: signed/unsigned "
  echo ""
}


while [[ $# -gt 0 ]]
do
    key="${1}"
    case ${key} in
    -o|--output)
	OUTPUT="${2}"
	shift
	shift
	;;
    -s|--sft)
        SFT="${2}"
        shift # past argument
        shift # past value
        ;;
    -n|--network)
        NET="${2}"
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
	echo -e "$(tput setaf 1)NO output path supplied!\n$(tput sgr 0)Please enter the directory of sampleTree.RData! (i.e. -o [path])"
	exit 1
fi

if [ -z "$SFT" ] || [[ ! "$SFT" =~ ^[0-9]+$ ]] ; then  # $SFT is number and not null
	echo -e "$(tput setaf 1)NO soft-threshold supplied!\n$(tput sgr 0)Please consult Part2_ScaleFree_topology.pdf! (i.e. -s [1~20])"
	exit 1
fi

if [ "$SFT" -gt 20 ] || [ "$SFT" -lt 1 ]; then
	echo -e "$(tput setaf 1)Chosen soft-threshod is out of range (1~20) $(tput sgr 0)"
	exit 1
fi

if [ "$NET" != "signed" ] && [ "$NET" != "unsigned" ] ; then
    $NET = "unsigned"
    echo "Network type is not supported by WGCNA, build unsigned network by default!"
fi

# run R
Rscript "${DIR}/WGCNA_part3.R" $OUTPUT $SFT $NET


