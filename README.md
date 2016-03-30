# README #

This analysis pipeline is designed to run a modified version of PyClone - allowing subclonal copy number states. This analysis relies entirely on PyClone. Therefore, this must be installed and set-up. 
An excellent set of instructions of how to get this up and running can be found at 
https://bitbucket.org/aroth85/pyclone/wiki/Home

This analysis pipeline was implemented in the paper https://www.ncbi.nlm.nih.gov/pubmed/?term=26940869. If you use this tool, please cite this paper. 

The analysis pipeline can work as a command line tool. The following input is required:

  <patient>              # The identifier for the patient/sample

  <saveDir>              # Where do you want to save the output?

  <ascatDir>             # Where is the ASCAT copy number data saved? (This should be a directory)

  <snvDir>               # What is the mutation table saved? (This should be a directory)

  <PyClone>              # Where is PyClone installed, e.g. "/farm/babs/redhat6/software/python/bin/PyClone"

  <template.config.yaml> # Where is the pyClone template? 

An example of each and every file can be found at https://bitbucket.org/nmcgranahan/clonalneoantigenanalysispipeline/downloads

An example implementation of the script is as follows:

Rscript ~/Documents/Work/clonalDissectionPyClone.R LMS025 ~/Documents/Work/AnalysisPipeline/ ~/Documents/Work/AnalysisPipeline/ASCAT/ ~/Documents/Work/AnalysisPipeline/SNV/ /farm/babs/redhat6/software/python/bin/PyClone ~/Documents/Work/AnalysisPipeline/template.config.yaml ~/Documents/Work/AnalysisPipeline/clonalDissectionFunctions.R