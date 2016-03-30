# README #

This analysis pipeline is designed to run a modified version of PyClone - allowing subclonal copy number states. This analysis relies entirely on PyClone. Therefore, this must be installed and set-up. 
An excellent set of instructions of how to get this up and running can be found at 
https://bitbucket.org/aroth85/pyclone/wiki/Home

The analysis pipeline can work as a command line tool. The following input is required:

  <patient>              # The identifier for the patient/sample

  <saveDir>              # Where do you want to save the output? (this should be a full path)

  <ascatDir>             # Where is the ASCAT copy number data saved? (a full path to a directory)

  <snvDir>               # What is the mutation table saved? (a full path to a directory)

  <PyClone>              # Where is PyClone installed,  e.g. "/farm/babs/redhat6/software/python-2.7.3/bin/PyClone"

  <template.config.yaml> # Where is the pyClone template? (a full path to the file)

An example of each and every file can be found at https://bitbucket.org/nmcgranahan/clonalneoantigenanalysispipeline/downloads

An example implementation of the script is as follows:

/usr/bin/Rscript ~/Downloads/ExampleFiles/clonalDissectionPyClone.R LMS025 ~/Documents/Work/AnalysisPipeline/ ~/Downloads/ExampleFiles/ASCAT/ ~/Downloads/ExampleFiles/SNV/ /farm/babs/redhat6/software/python/bin/PyClone ~/Downloads/ExampleFiles/template.config.yaml ~/Downloads/ExampleFiles/clonalDissectionFunctions.R

This analysis pipeline was implemented in the paper https://www.ncbi.nlm.nih.gov/pubmed/?term=26940869. If you use this tool, please cite this paper. The results on each tumour used, can be found at: https://bitbucket.org/nmcgranahan/clonalneoantigenanalysispipeline/downloads