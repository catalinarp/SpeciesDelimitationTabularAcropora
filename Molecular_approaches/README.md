# <b>Molecular species delimitation</b>

This folder contains the files and commands used for molecular species delimitation both using target-capture sequencing data and Sanger sequencing of PCR-based amplificated markers derived from target enrichment.

## Target-capture sequencing
Target enrichment data for the <b><i>Screening of target-enriched loci</i></b> phase, including contig fasta files for each sample that matched the baits designed both for the exons and the UCEs sets.

## Sanger sequencing
Chromatograms and alignments used for molecular analyses for stages: <b><i>Preliminary screening of available molecular markers</i></b> and <b><i>Implementation of target-enrichment derived markers in molecular species delimitation</i></b> 

## Downstream analyses
Folders containing scripts and examples of commands used for the analyses performed in each one of the stages detailed in the following table: 
<br>

<b>Summary Table | Techniques, loci and methods used in the different stages of the molecular analyses performed in this study.</b> Detailed information about the techniques, the number of loci, the number of individual samples, and the general pre-processing steps and downstream analyses used in each stage of the molecular approaches used in this study are displayed.<br>

| Stage                                                                                 | Molecular technique                                                                                                              | No. loci/markers [n= samples]                     | Pre-processing                                                                                                                           | Downstream analyses                                                                                                                               |
|---------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| <b><i>1. Preliminary screening of available molecular markers</i></b>                                  | PCR-based amplification followed by Sanger sequencing                                                                            | Three genetic markers (AcroCR, PMCA, FZD) [n= 36] | Chromatograms edition, sequence alignment and phasing                                                                                    | Genetic clustering, genetic distances and gene trees                                                                                              |
| <b><i>2. Screening of target-enriched loci</i></b>                                                     | Target enrichment and high-throughput sequencing of conserved elements (exons and UCEs) captured using the hexacoral-v2 bait set | 2060 loci (1026 exons, 1034 UCEs) [n= 9]          | Reads de-multiplexing and trimming, contigs assembly and probe matching followed by <b>1)</b> Phasing Laninsky pipeline or <b>2)</b> Phasing PHYLUCE pipeline  | For <b>1)</b> Genetic clustering (1889 loci), SNAPP species tree (210 loci); and <b>2)</b> Allele sharing-based approaches and extended species trees (79 loci) |
| <b><i>3. Implementation of target-enrichment derived markers in molecular species delimitation</i></b> | PCR-based amplification followed by Sanger sequencing                                                                            | Three genetic markers (TDH, DOPR, ASNA) [n= 36]   | Chromatograms edition, sequence alignment and phasing                                                                                    | Genetic clustering, genetic distances, gene trees, species trees, coalescent and allele sharing-based approaches                                  |
