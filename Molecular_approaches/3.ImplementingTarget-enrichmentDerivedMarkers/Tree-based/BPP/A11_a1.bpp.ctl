          seed =  -1

       seqfile = ~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/BPP_analyses/3TC.phy
      Imapfile = ~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/BPP_analyses/3TCmap.txt
       outfile = 3TC_run2.out
      mcmcfile = 3TCmcmc_run2.txt
      
  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2   * species delimitation rjMCMC algorithm0 and finetune(e)
 speciesdelimitation = 1 1 2 1  * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  A  B  C  D
                    12  12 12  1
                 (((A, B), C), D);
        diploid =   1  1  1  1  * 0: phased sequences, 1: unphased diploid sequences
                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 3  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.004   # invgamma(a, b) for theta
      tauprior = 3 0.002    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 20000
      sampfreq = 1
       nsample = 200000
