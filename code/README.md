# README.md

## Code

### Notebooks

Notebooks are used to either test tools or perform analysis.

- `Motif_Scoring_and_Extraction.ipynb`: Larger notebook that acts as a place to manually try out all the fucntions that go into running our pipeline.
- `run_all_alignments.ipynb`: Python notebook used to run MSE.pipeline() on all 3,542 alignment files. Successfully runs all alignment files through the pipeline and writes each output to a CSV. Currently only outputting the "filtered" table portion of the output.
- `mapping_TFBS_10Sept2018.ipynb`: testing if the `pipeline()` function can be used in a loop. This will need to be run across all the alignment files.
- `zelda_20Aug2018.ipynb`: This is attempting to run the tools on isolated files.
- `zelda_12July2018.ipynb`: This is attempting to run the tools on isolated files. 

### Utilities / Tools

These are the files that house the functions that are used in the notebooks. They are largely built to make the `pipeline()` function work. 

- `threshold_setter.py`: One of the main problems with consistantly mapping TFBS is that in order to say a TFBS is present, there need to be a score threshold set. We formed our own criteria that sets this threshold.
- `MSE.py`: These are the functions that aid in the actual mapping of TFBS once the threshold has been set. Contains the main `pipeline()` function that takes in an alignment file, a motif, and a threshold, and outputs tables containing the TFBS's and corresponding information.