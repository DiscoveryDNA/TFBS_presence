# README.md

## About 

The purpose of this project is to see if Transcription Factor Binding Sites (TFBS) presence is correlated with enhancer function. This repository function is to build a reproducible way to map motifs onto TFBS across alignemnts. In order to achieve this we need a realiable pipeline that gives a score to each position of an alignment file. 

The program essentially:

1. Reads in an alignment 
2. Splits each sequence in the alignment into its own biopython sequence object with removed gaps (dictionary) and creates a dicitionary reference that allows each seperated raw sequence postion to be mapped back to the original alignment position. 
3. Uses Biopython's motif package, the program gives a score on each raw sequence position corresponding to how likely an input Position Weight Matrix (PWM) is found on that sequence. This is run on both the positive and negative strand.
4. The output is several tables mapping each score, at each position, across all sequences, with the raw and alignment postion. 

In order to set a score threshold for determing presence and absence at each positition:

1. A randomly generated nucleotide fasta file is created and the distribution of values is decided.
2. ect


## Code

### Notebooks

Notebooks are used to either test tools or perform analysis. 

- **Current Work** `mapping_TFBS_10Sept2018.ipynb`: testing if the `pipeline()` function can be used in a loop. This will be need to run across all the alignment files.
- `Motif_Scoring_and_Extraction.ipynb`: Larger notebook that acts as a place to manually try out all the fucntions that go into running our pipeline. Functions in this notebook may be outdated from what is used in the actual tools run. 
- `zelda_12July2018.ipynb`: This is attempting to run the tools on isolated files. 
- `zelda_20Aug2018.ipynb`: This is attempting to run the tools on isolated files.

### Utilities / Tools

These are the files that house the furtnions that are used in the notebooks. They are largely built to make the `pipeline()` function work. 

-`threshold_setter.py`: One of the main problems with consistantly mapping TFBS is that in order to say a TFBS is present, there need to be a score threshold set. We formed out own criteria that sets 
- `MSE.py`: These are the functions that aid in the actual mapping of TFBS on 


## Workflow 

Below is a schematic of the entire program workflow. 

![workflow](https://github.com/DiscoveryDNA/TFBS_presence/blob/master/img/TFBS-workflow-01.jpg)


## Resources

[Sync Upstream](https://help.github.com/articles/syncing-a-fork/)

