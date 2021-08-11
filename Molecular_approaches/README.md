# <b>Molecular species delimitation</b>

This folder contains the files and commands used for molecular species delimitation both using target-capture sequencing data and Sanger sequencing of PCR-based amplificated markers derived from target enrichment.

## + Target-capture sequencing


## + Sanger sequencing
Chromatograms and alignments used for molecular analyses for stages: <b>Preliminary screening of available molecular markers</b> and <b>Implementation of target-enrichment derived markers in molecular species delimitation</b> (see following table).

<br>
<b>Summary Table | Techniques, loci and methods used in the different stages of the molecular analyses performed in this study.</b> Detailed information about the techniques, the number of loci, the number of individual samples, and the general pre-processing steps and downstream analyses used in each stage of the molecular approaches used in this study can be found in the following table:
<br>

| Stage                                                                                 | Molecular technique                                                                                                              | No. loci/markers [n= samples]                     | Pre-processing                                                                                                                           | Downstream analyses                                                                                                                               |
|---------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| <b>Preliminary screening of available molecular markers</b>                                  | PCR-based amplification followed by Sanger sequencing                                                                            | Three genetic markers (AcroCR, PMCA, FZD) [n= 36] | Chromatograms edition, sequence alignment and phasing                                                                                    | Genetic clustering, genetic distances and gene trees                                                                                              |
| <b>Screening of target-enriched loci</b>                                                     | Target enrichment and high-throughput sequencing of conserved elements (exons and UCEs) captured using the hexacoral-v2 bait set | 2060 loci (1026 exons, 1034 UCEs) [n= 9]          | Reads de-multiplexing and trimming, contigs assembly and probe matching followed by 1) Phasing Laninsky pipeline or 2) Phasing PHYLUCE pipeline  | For 1) Genetic clustering (1889 loci), SNAPP species tree (210 loci); and 2) Allele sharing-based approaches and extended species trees (80 loci) |
| <b>Implementation of target-enrichment derived markers in molecular species delimitation</b> | PCR-based amplification followed by Sanger sequencing                                                                            | Three genetic markers (TDH, DOPR, ASNA) [n= 36]   | Chromatograms edition, sequence alignment and phasing                                                                                    | Genetic clustering, genetic distances, gene trees, species trees, coalescent and allele sharing-based approaches                                  |
