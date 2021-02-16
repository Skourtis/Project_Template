# Project_Template
Template

Hi Maria, 


This is a Project Template I have created for Git.
To get started using Rstudio, create a new project, and select the option to clone something from Git. Name the project anything your heat desires. 
Give it the address of the is Repository on Git and once everything is downloaded, Run script 01 (the script should be open in the the project you just created).
Running the script, should install all the packages and download all the files that we will be using (and so it's more reproducible as we are using the same versions of tools).

If you made it here without any major errors and problems amazing!

There is a also a 2020MQ044.Rmd file. This will be our main working space and output. 
Trying getting yourself familiar with the concept of Rmarkdown (which is similar to Jupiter Notebooks in Python)

I have pre-written some functions in functions.R and they are being called by the markdown. 
Unzip the Processed Proteomics data in Datasets -> Raw. Unzip this.
The Rmarkdown will look for the unzipped folder and you can then produce the PTXQC quality control report as a first step. 
Then using the following lines of code you should get matrix of Samples against protein Uniprot IDs with protein Abundances.

This is where the Data analysis fun starts! 
I hope you've made it to here. (somethings things just don't work but it's nice when they do).

So for week 1 lets aim to:
	1. Set up our Git communication
	2. Produce the final QC
	3. Get a general idea of how the samples are behaving
	4. impute (remove proteins which are rare) and find Differentially expressed Genes
	5. Cluster the proteins and samples. 
	6. Do some simple enrichment analysis on those. 