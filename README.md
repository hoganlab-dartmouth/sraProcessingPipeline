# sraProcessingPipeline
sra Processing Pipeline

Python scripts in LocalScripts are test scripts that were written to run locally. Python scripts in the DiscoveryScripts were written to run on Discovery. Python scripts in the Graphing folder are for graphing.

# steps to run pipeline
1. Clone git repo to discovery home folder
2. Go into the scripts and change out f004ky5 to your username, you may need to change more than just this to make sure scripts are user specific. I would do a search for f00 and then change out anywhere that string appears.
3. qsub indexer.pbs
4. conda load anaconda3
5. python run_lite.py
6. qsub csv_builder.pbs
7. output csvs should appear in git repo folder
