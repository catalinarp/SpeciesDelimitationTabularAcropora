# <b>Sanger sequencing data</b>

Sequences were obtained by PCR-based amplification using both forward and reverse primers from 36 tabular Acropora samples collected in Sesoko Island (Okinawa, Japan). The markers used were:

#### Previously reported markers:
+ Mitochondrial putative control region (AcroCR) flanked by 12S ribosomal RNA at the 5' end and the gene cox3 at the 3' end
+ Plasma membrane calcium-transporting ATPase (PMCA) as published by Ladner & Palumbi, <i>Mol. Ecol.</i> 2012, 21(9):2224-38.
+ Frizzled-4 like homolog (FZD) extended from the published by Ladner & Palumbi, <i>Mol. Ecol.</i> 2012, 21(9):2224-38.

#### Markers defined from loci captured using target enrichment strategies:
+ L-threonine 3-dehydrogenase (TDH)
+ Dopamine receptor 2-like (DOPR)
+ ATPase ASNA1 homolog (ASNA)

<b>Primers and conditions for PCR-based amplification and Sanger sequencing.</b> General PCR conditions: start 1 sec 95°C; 1 min 95°C; [30 sec 95°C; 30 sec Tº annealing (Ta); 2 min 72°C]x Number of cycles; 10 min 10°C. *Difficult samples were re-amplified and sequenced using M13-tailed PCR primers and M13 tails, respectively. M13-tail-F: TGTAAAACGACGGCCAGT, M13-tail-R: CAGGAAACAGCTATGAC. a) Product length and primers from previous studies (Ladner & Palumbi, Mol. Ecol. 2012, 21(9):2224-38).b) Product length obtained from primers designed in-house. c) Product length was extended by complementing sequences obtained using previously reported primers (Ladner & Palumbi, Mol. Ecol. 2012, 21(9):2224-38), with sequences obtained using primers designed in-house. GenBank accession numbers for the allele sequences obtained with each marker are also shown (GenBank IDs).

|            Loci (GenBank IDs)              |                                                       PCR primers (5’ - 3’)                                                      |          PCR conditions          |                                           Sequencing primer     (5’ - 3’)                                        |      Product length (bp)      |
|:------------------------------------------:|:--------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------:|:----------------------------------------------------------------------------------------------------------------:|:-----------------------------:|
|     AcroCR     (MT945838 - MT945873)       |     F:   GCCCCTCAAGAGGGTTTCTA     R:   CTAGACAGGGCCAAGGAGAAG                                                                     |     Ta: 55º     55   cycles      |     F: PCR primer F     R: PCR primer R                                                                          |     1265   – 1352b            |
|     PMCA     (MT945609 - MT945656)         |     Fa:   AAGGAATTGGTGGCTTTCCT      Ra:   CACAGACGACCATCTTTCCA                                                                   |     Ta:   53º     50   cycles    |     F:   GAATTGGTGGCTTTCCTGAG     R:   CGACCATCTTTCCACTACCTTC                                                    |     545a                      |
|     FZD     (MT945657 - MT945718)          |     F1a:   TATGGCTGCGACAATTTGGT     R1a:   GCTAGCGTTTCGAGTTCCAC     F2:   CCTTGAGTTGGTTCCTTGCT     R2:   CGCCTAGACAGCAGCTAAAA    |     Ta:   55º     50   cycles    |     F1:   CCTTGAGTTGGTTCCTTGCT     R1:   TCGAGTTCCACCGTTCTTCT     F2:   PCR primer F*     R2:   PCR primer R*    |     639a     994   – 1006c    |
|     TDH     (MT945719 - MT945777)          |     F:   TTTTTCTTTCACTTTTGGCTGT     R:   ATCTCTGCTGCAATCCCAAT                                                                    |     Ta: 53º     50   cycles      |     F: PCR primer F*     R: PCR primer R*                                                                        |     736   – 744b              |
|     DOPR     (MT945778 - MT945837)         |     F:   AGGGTCAGGTTTTTGGGAAT     R:   GAGTTTTGACCGTCAGTTGG                                                                      |     Ta: 53º     50   cycles      |     F: PCR primer F*     R: PCR primer R*                                                                        |     747   – 760b              |
|     ASNA     (MT945874 - MT945940)         |     F:   CTGTGTGCTGGCGAAAAA     R:   GAAAGGCCCCTCTATTTTCA                                                                        |     Ta: 53º     50   cycles      |     F: PCR primer F*     R: PCR primer R*                                                                        |     748   – 763b              |

#### Files
+ Sequencher v4.0 files (.SPF) containing chromatograms and contigs assembled from Sanger sequencing data (.ABI files)
+ Mafft v7.471 alignments of consensus contigs (consensus.mafft.fasta) and of phased haplotypes (mafft.alleles.fasta).

Forward and reverse chromatograms of each individual were assembled into contigs in Sequencher. For each homozygote, the contig was checked and cleaned then received the name of the individual sequenced (18OkiXX). For length-variant heterozygotes, the double peaks were called either manually or using Sequencher's "Call Secondary Peak" function, then the two haplotypes were reconstructed using the program Champuru v1.0 (http://jfflot.mnhn.fr/champuru/). The contig received the name of the individual (18OkiXX-CHPR), and the reconstructed haplotypes of each individual received the name of this individual followed with underscore "a" and "b", respectively. For heterozygotes having two alleles of equal lengths, we used SeqPHASE (https://eeg-ebe.github.io/SeqPHASE/index.html) and PHASE to infer their haplotypes. The contig of each such individual received the name of this individual (18OkiXX), and its two inferred haplotypes received its name followed with underscore "a" or "b", respectively.

