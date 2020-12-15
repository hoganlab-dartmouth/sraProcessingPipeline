# sraProcessingPipeline
sra Processing Pipeline

Python scripts in LocalScripts are test scripts that were written to run locally. Python scripts in the DiscoveryScripts were written to run on Discovery. Python scripts in the Graphing folder are for graphing.

# steps to run pipeline
1. Clone git repo to discovery home folder
2. Go into the scripts and change out f004ky5 to your username, you may need to change more than just this to make sure scripts are user specific. Fastes way to do this would be to do a search for f00 and then change out anywhere that string appears. I think this needs to be done for every file (including pbs scripts). 
3. qsub indexer.pbs
4. conda activate anaconda3
5. python run_lite.py
6. qsub csv_builder.pbs
7. output csvs should appear in git repo folder

# location of Salmon calls and fastq-dump call for changing parameters. 
1. Salmon indexer call appears in indexer.py line 33. Using the default parameters here, but should change to reflect average read length of sample set. Current setting is optimized for a mean of 75bp.
2. SRAToolkit fastq-dump call appears once in quant.py line 78. Parameters: --gzip, --skip-technical, --readids, --split-e (this is the same as split-3 I think), --clip. 
https://edwards.sdsu.edu/research/fastq-dump/
3. Salmon quant call appears twice in quant.py line 90 and 97. These should be combined into a single. Make sure if you change one to change the other. The only parameters we've changed here is -l A for library type: Automatic and -r input_file for unpaired inputs. Since we're doing independent paired reads, we should specify --fldMean and --fldSD once we know the standard deviation (gaussian) and mean of our sample set. 
