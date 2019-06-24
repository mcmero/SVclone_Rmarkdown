# SVclone_Rmarkdown

Generate plots and tables from the SVclone paper and supplementary material.

### How do I get set up? ##

Install [R](https://www.r-project.org/), [Rstudio](https://rstudio.com/products/RStudio/) and the following packages:

    install.packages('grid')
    install.packages('gridExtra')
    install.packages('data.table')
    install.packages('RColorBrewer')
    install.packages('ggplot2')
    install.packages('circlize')
    install.packages('survival'
    install.packages('ggfortify')
    install.packages('pROC')
    
    source('https://bioconductor.org/biocLite.R')
    biocLite('GenomicRanges')
    biocLite('biomaRt')

Into the same directory as this repository, clone the PCAWG colour palette code:
    
    git clone https://github.com/ICGC-TCGA-PanCancer/pcawg-colour-palette.git
    
The simulation_results.Rmd notebook can now be generated. The simulation data can also be regenerated by running the [SV simu pipeline](https://github.com/mcmero/sv_simu_pipe).

To generate the metrics.Rmd notebook, the 001 metastases from doi:10.1038/ncomms7605 must be obtained. Instructions for subsampling and creating the in-silico sample mixtures must be followed as described in the SVclone paper (https://doi.org/10.1101/172486). The script for generating the _in silico_ mixtures can be found in this repository (`make_insilico_mixtures.sh`) and requires [samtools](https://sourceforge.net/projects/samtools/files/samtools/). The following files are required for the original bM and gM samples:

* variant calls in Mutect's call stats format (put under 001_truth/001*M_SS_mut_merged.call.stats.txt)
* Battenberg output (subclones files under 001_truth/001*M_subclones.txt)
* svinfo output from SVclone's count step (under 001_truth/001*M_merged_list_svinfo.txt)
    
For each mixture, bases must be counted at each variant locus, and added to the 001_snv_counts folder, with the following name convention (001bM_p\<p1\>_001gM_p\<p2\>_merged_sorted_reheader.basecounts.vcf), where p1 and p2 are the proportions (e.g. 001bM_p0.9_001gM_p0.1_merged_sorted_reheader.basecounts.vcf). SVclone (all steps) and Pyclone must be run on all mixtures, with the output for SVclone under 001_mix_results/001bM_p\<p1\>_001gM_p\<p2\> and Pyclone 001_pyclone_out/001bM_p\<p1\>_001gM_p\<p2\>. Battenberg must be run on each mixture and the corresponding output (subclones file) should be added to a `battenberg_out` directory at the base of the repository. 

To generate the notebook containin PCAWG analyses, access to the ICGC/TCGA PCAWG is required. The following files are required (see `analyse_select_pcawg.Rmd` for file names):

* PCAWG sample sheet
* expanded release spreadsheet
* annotated purities and ploidies
* clinical donor spreadsheet (PCAWG)
* ICGC clinical donor information
* PCAWG amplified FBI fraction
* specimen histology spreadsheet
* list of PCAWG driver genes from previous knowledge
* list of PCAWG SNV coding hits in driver genes

Obtain SNV calls, as well as both full and consensus CNAs, checking that the system paths correspond to the notebook paths. Obtain SVclone results from the PCAWG jamboree. The notebooks can now be run.
