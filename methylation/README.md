Based on prior QC work by Christian M. Page and Haakon E. Nustad and
the bhklab/IlluminaEPICmethylation pipeline on github by Christopher Eeles and Joanna Pryzbyl

set up steps:
steps that need to be done, specific to this pipeline on TSD:
- navigate to "/tsd/p581/data/durable/Projects/MoBa_genomics/moba_qc_pipeline"

- run:
$ source .venv/bin/activate
This activates the virtual environment, enabling snakemake

if you wish to jump out of the virtual enviroment, do:
$ deactivate

when R is run, the packages specific to this project are loaded. These are listed in the renv.lock file. If another package must be installed, install.packages("X", repos = "repository_path") in an interactive session should install the package from the correct resources on TSD (the cran and bioconductor copies). the repository path must be specified, and the available options can be seen in the top of the renv.lock file. One can then update the renv.lock file from inside R by using renv::snapshot(). If one interactively install packages, and those packages creates troubles, do a renv::restore() to reset and reload packages from renv.lock.

adding new steps:
When adding scripts to the pipeline, look at the structure of an existing R-script. An important step is the source(".Rprofile"), which enables the local renv-library. Similarly, look at the Snakefile how the rules are defined.

Pre-quality control:
important pre-step: the sample sheet needs the beadchip array column to have column name: Sentrix_ID, 
and the position on the array to have column name: Sentrix_Position
This is important for minfi to automatically recognize correct columns to identify the .idat files.

All samples you wish to run QC for must be listed in your sample_sheet.csv, which must be located in the folder you specify in config.yaml. It should be placed together with the .idat files. This .csv file is recognized automatically with the minfi function applied to read in the sample_sheet. A current step implemented detects if you have a mix of children and older individuals, which will create an error with an error message, encouraging the user to split his/hers samples into two .csv files, and running these seperately. This is because the quality of the QC will increase if such biologically related differences are not present.

If your data is larger than 500 EPIC samples or 850 450k samples, the current TSD node with 60GB RAM will run into memory issues. This also requires a split of your data.

 


Quality Control Pipeline:

summary:
- load iDat files to RGset (all information in minfi object format)
- check spread for all samples in PCA plot of approx. 600 control probes (plot is created for manual inspection. Look for outlying samples. compares information across samples)
	- possibly create plot for some specific control probes and/or median intensity across control probes (plot is created for manual inspection. Look for outlying samples)
- drop samples with less than 90% (or possibly some other cutoff) probes detected. (detection p value larger than 0.01 for more than 10% of probes -> remove sample)
- drop samples with bisulphite conversion rate less than 90% (this cutoff can be set down to 80%, but should not be set any lower)
	This is calculated with bscon(), which uses bisulphite conversion control probes to estimate conversion rate. algorithm calculates within each sample
- single-sample-Noob method used for background subtraction and dye-bias correction. This algorithm also only uses information within each sample to correct each sample. Has been shown to be one of the better background subtraction and dye-bias correction methods. 
-  create some plots comparing methylation values before and after ssNoob.
- filter out poor quality probes; probes shown to have a high detection p value (higher than 0.01 (default), 0.02 or 0.05 the usual cutoffs) for more than 0.05 (cutoff can be changed) of the samples. naturally, this step compares information across samples.
- remove probes shown to be cross-hybridizing, non-CpG methylation (CpH) and probes influenced by SNPs close to interrogation CpG site. These can be removed with a blacklistlist of known CpG sites.
	cross-hybridising: probes shown to possibly have more than 1 binding site in the genome.
	Some probes designed to interrogate methylation at CpH (H = A, C or T) sites.
	SNP influenced: probes affected by a SNP at or close to the target site. usual nearness cut-off: 3 or 5 base pairs from interrogation site. These CpGs will be highly influenced by SNP variant, both in terms of methylation value and chemical binding to probes.
- BMIQ: probe type normalisation. the Beadchips have 2 different probe designs, I and II. These have a clear difference in dynamic range for the obtained methylation measurements. BMIQ projects the global distribution of probe type II measurements onto the distribution of probe type I, meaning that the methylation values from probe type II measuements are slightly shifted. This is especially a usual step when analysing regions of data, i.e. looking for differentially methylation regions. It is also believed that the dynamic range of probe type I is more correct than probe type II. BMIQ also calculates only within each sample.

- the resulting data ready for delivery.





