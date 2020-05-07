## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Usage](#usage)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## About The Project

In this project, I use WGCNA to mine the gene expression data extract from RNA-seq dataset of Leucegene corhort. Objectif is to discover new gene expression signature or biomarkers related to AML.

* [WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)

Weighted gene co-expression network analysis (WGCNA), is a widely used data mining method especially for studying biological networks based on pairwise correlations between variables. It has been most widely used in genomic applications (While it can be applied to most high-dimensional data sets). 
WGCNA allows to define modules (clusters), intramodular hubs, and network nodes with regard to module membership, to study the relationships between co-expression modules, and to compare the network topology of different networks (differential network analysis). WGCNA can be used as a data reduction technique, as a clustering method (fuzzy clustering), as a feature selection method, as a framework for integrating complementary (genomic) data (based on weighted correlations between quantitative variables), and as a data exploratory technique.

* [Leucegene](https://leucegene.ca/)

Quebec-based team of multidisciplinary scientists, clinicians, economists and lawyers with the mission to change the treatment of acute myeloid leukemia (AML) using precision therapy approaches. The AML cohort provides by far one of the largest genomic datasets for AML studies.


### Built With
* [R 3.6.1](https://www.r-project.org/)
* [WGCNA package](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)
* [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)
* [KEGGREST](https://bioconductor.org/packages/release/bioc/html/KEGGREST.html)



<!-- USAGE EXAMPLES -->
## Usage
### use in Rstudio
Please run the R code by order, each part will generate informatif plots for user to make parameter decisions for the following part. Each part will automaticlly save a Rdata for the following part.

* WGCNA_1_outlierDetection.R : 
  remove low variance genes and genes with too much missing values, detect outlier samples

* WGCNA_2_scaleFreeTopology.R : 
  Gather topology information for soft-threshold selection

* WGCNA_moduleDetection.R :
  Detect gene expression modules, calculate correlation between modules and traits
  
* WGCNA_moduleAnalysis.R :
  Analyze module of interest, calculate gene-trait correlation, find hubgene, GO enrichment analysis, KEGG path way analysis

### use in Terminal
Please run the bash code by order, each part will launch the corresponding R code. Each part will automaticlly save a Rdata and generate informative pdf plots for parameter selection in the following part.

* WGCNA_part1.sh:  
 remove low variance genes and genes with too much missing values, detect outlier samples  
 ***Usage***: ./WGCNA_part1.sh  **-o**  [Output Path]  **-d**  [Expression data]  **-t**  [Traits data]  **-l**  [log base] **-c**  [add constant before log]
 ***Help***:  ./WGCNA_part1.sh  **-h**
 
* WGCNA_part2.sh:  
  Gather topology information for soft-threshold selection  
  ***Usage***:   ./WGCNA_part2.sh  **-o**  [PATH of sampleTree.Rdata]  **-c**  [Height]  
  ***Help***:    ./WGCNA_part2.sh  **-h**

* WGCNA_part3.sh:  
  Detect gene expression modules, calculate correlation between modules and traits  
  ***Usage***: ./WGCNA_part3.sh **-o** [PATH of topology.Rdata] **-s** [Soft-threshold]  
  ***Help***:  ./WGCNA_part3.sh **-h**

* WGCNA_part4.sh:  
  Analyze module of interest, calculate gene-trait correlation, find hubgene, GO enrichment analysis, KEGG path way analysis  
  ***Usage***: ./WGCNA_part4.sh **-o** [PATH of network.Rdata] **-m** [Module] **-t** [Trait] **-i** [GeneID type]  
  ***Help***:  ./WGCNA_part4.sh **-h**




<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

[Ryan Xiao Ju](https://twitter.com/RyanXJu0505) - ryan.ju0505@gmail.com

Project Link: [https://github.com/RyanXJu/WGCNA/](https://github.com/RyanXJu/WGCNA/)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [IRIC Bioinformatics platform](https://www.iric.ca/en/research/platforms-and-infrastructures/bioinformatics)
* [Leucegene](https://leucegene.ca/)




