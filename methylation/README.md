<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Introduction](#introduction)
- [Setup](#setup)
    - [TSD](#tsd)
        - [Start interactive R:](#start-interactive-r)
        - [exit R](#exit-r)
- [Configuration](#configuration)
    - [Global config file globalConfig.yaml](#global-config-file-globalconfigyaml)
    - [Local config-files](#local-config-files)
- [Pre-quality control:](#pre-quality-control)
- [Running the pipeline](#running-the-pipeline)
    - [Removing bad samples](#removing-bad-samples)
    - [Finalizing the run](#finalizing-the-run)
    - [Canned wrapper script](#canned-wrapper-script)
- [Results](#results)
    - [Final](#final)
    - [Logs](#logs)
    - [Plots](#plots)
    - [Qc_results](#qc_results)
- [Quality Control Pipeline documentation/descrition](#quality-control-pipeline-documentationdescrition)

<!-- markdown-toc end -->

# Introduction
Based on prior QC work by Christian M. Page and Haakon E. Nustad and
the bhklab/IlluminaEPICmethylation pipeline on github by Christopher Eeles and Joanna Pryzbyl.
Developed by Hanna Zdanowicz, Haakon E. Nustad and Gutorm Thomas Høgåsen.

# Setup
General set-up steps (outside of TSD): Clone the git repository. Every
command must be run from inside the main directory.

Then install (if not installed)
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
An exported environment (suitable for running the pipeling) is available
on Setup/metQc-env.yml or metQc-env.txt .  Among many things it will
install R, python and snakemake. The txt-version was created by the
command\ conda list --explicit > metQc-env.txt\

It is suitable for a linux-64 platform and can be used to create an
environment (provided conda is installed) by\
conda create --name metQc --file metQc-env.txt


The yml version on the other hand is more generic and made by
conda env export > metQc-env.yml\

It should work to create a working environment with 

conda env create --file metQc-env.yml

Since you use conda, you do not want to use renv (see the TSD section below),
so make sure renv/, renv.lock and .Rprofile are located in Setup
folder and make sure that in input/globalConfig.yaml the renv variable
is set to to FALSE.

## TSD
Set up on TSD, needs a litlle extra tweaking. 

Unpack the zipped pipeline folder at preferred destination on TSD, move within the folder and use python3 to set up virtual environment and install snakemake:
$python -m venv .venv
$pip install snakemake

Move renv/, renv.lock and .Rprofile out of Setup folder to root folder.

Change globalConfig.yaml renv variable to TRUE.

For installation of packages into local R library:
Start interative session with R from root pipeline folder. This will trigger renv, which is version control and package management software for R.
Run the following commands:

### activate the virtual env:
$ source .venv/bin/activate
This activates the virtual environment, enabling snakemake.

(If you wish to jump out of the virtual enviroment, do:
$ deactivate
)

### Start interactive R:
R

#trigger dependency installation if initialization did not:
renv::restore()

### exit R
q()


# Configuration

## Global config file: globalConfig.yaml
Edit the globale config file globalConfig.yaml, it must be customized
for your environemnt. 

Important variables:
- root_data_path points to where idat-files and samplesheets are
  located. The exact locations are found in the local config-files.
- ouput_path is where all the results logs and whatnot will end
  up. The pipeline might create terrabytes of data depending on your
  input
- max_sample_size_* Unless you have loads of memory, you want to
  reduce the size of samples being processed in one go. See the
  documentation and comments in the yaml-file for more details.
- renv If you cannot use conda as described below, an alternative
  using pip & renv is possible.

## Local config-files
Each dataset has one or more config-file (called a local config file)
that can be used to override the global one.

If you use the standard/original 'MoBa' path-structure you don't need
to change anything.

Go through the local config.yaml file (usual names met003.yaml,
met004.yaml, etc) in input folder and change parameters to your
specific file paths and names, especially data directory (idat and
sample sheet file path), data type (450k or EPIC). 

Additionally, add file in input folder:
manual_remove_samples_bad_density.txt. This file is meant to contain
sample ids of samples you want to remove based on having outlying
genomewide DNA methylation patterns. This should be filled with sample
ids based on output from the first and second part of the
pipeline. Good to go!


# Pre-quality control:
Important pre-step: the sample sheet needs the beadchip array column to have column name: Sentrix_ID, 
and the position on the array to have column name: Sentrix_Position.
This is important for minfi to automatically recognize correct columns to identify the .idat files.

All samples you wish to run QC for must be listed in your
sample_sheet.csv, which must be located in the folder you specify in
config.yaml. It should be placed together with the .idat files. This
.csv file is recognized automatically with the minfi function applied
to read in the sample_sheet. A current step implemented detects if you
have a mix of children and older individuals, which will create an
error with an error message, encouraging the user to split their
samples into two .csv files, and running these seperately. This is
because the quality of the QC will increase if such biologically
related differences are not present.

# Running the pipeline

Most users will want to use the [Canned wrapper
script](#canned-wrapper-script) described below. The sections below
describe the different partst of the pipleline and how to run them with snakemake commands.

## Removing bad samples
The pipeline is organized in 3 parts, where one should inspect the
plots produced located in plots/ after the first and second part.

Based on the plots produced in first and second part, you should add sample IDs to the file input/(x)_manual_remove_samples_bad_density.txt if you think they should be removed. Remove samples showing especially outlying patterns compared to the rest, based on the boxplot of control_probes and genome wide density plot.

The pipeline can be run in several ways, with wrapper shell scripts or direct snakemake commands, depending on the users preferences. But in general, there is 1 command that should be run:
$ snakemake --core 1 --configfile input/globalConfig.yaml input/met001.yaml --snakefile Snakefile third_part

This command runs the entire pipeline if the input/(x)_manual_remove_samples_bad_density.txt file exists. If it is the first time the pipeline is run, one should first run the first_part and second part, then inspect the plots to consider putting some sample IDs in the input/(x)_manual_remove_samples_bad_density.txt file. This is done with the following commands:

$ snakemake --core 1 --configfile input/globalConfig.yaml input/met001.yaml --snakefile Snakefile first_part
$ snakemake --core 1 --configfile input/globalConfig.yaml input/met001.yaml --snakefile Snakefile second_part

## Finalizing the run

After one has run the third_part, the following command will clean up the output files, removing unnecessary files:
$ snakemake --core 1 --configfile input/globalConfig.yaml input/met001.yaml --snakefile Snakefile clean_up

## Canned wrapper script 
The wrapper script runOneorMore.sh will make the process of running
the pipeline easier, especially if your resources are limited. 
In
globalConfig.yaml, you need to specify 
- how many samples you can run in
memory at once. 2 estimates are given to make this easier. If your
amount of samples are larger than this amount, the pipeline will
return an error message the first time it runs, but will create
temporary config files where the data is split into batches.
- what part of the pipeline you want to run - edit the rule parameter.
This should in most cases be third_part.

This process is entirely controlled if you use the wrapper
script, but you will have to edit it describing what step you are running.

If you wish to run one local config file, met001.yaml for example, you
do this with the following command:

bash runOneorMore.sh input/met001.yaml

you can give several input config files if desired (the syntax assumes
bash):

bash runOneorMore.sh input/met00[147].yaml 

If the amount of samples in a set are more than the max_sample_size
specified in globalConfig.yaml, the data will be split automatically
into met001-1, met001-2, met001-3, etc. and run for each if these
batches. The scripts can be modified to allow more cores and which
rule to run. This can be changed by changing the variables inside the
scripts: cores rule

No arguments should be given to runOneorMore.sh if one wants to run through all local config files listed in the input folder. It is organized to run for all config files with the name consisting of "met" and ".yaml". Therefore, if the input folder has met001.yaml, met002.yaml, met003.yaml, it wil run for each of these. If any of these contains more samples than your hardware allows (max_sample_size), it will create subbatches. This script can also be modified with cores and rule. 

# Results
All results are stored within the Runs folder, with the analysis_name given as a subdirectory. The analysis_name is the name of the config file, minus the .yaml part. This is set automatically. 

## Final 
The end matrices with data:

This section has been moved completely to https://github.com/folkehelseinstituttet/mobagen/wiki/Methylation#QC

## Logs
Additional files with valueable information:
logs/:
- The Rout files are log-files containing information about the R session and the functions applied to the data. They contain information of what is done to the data during the different scripts. Additionally, if an error occurs in a script, the error message is given in the associated log file. These files are primarily meant for debugging of the pipeline, and can be inspected to understand the order of when the different functions and QC steps have been applied. 

## Plots
plots/:

{analysis_name}_3_boxplot_control_probes.pdf (boxplot of log2()-values from intensity measures from control probes for each sample)
{analysis_name}_3_control_probe_PCA_plot.pdf (PC2 against PC1 from PCA of log2()-values from intensity measures from control probe)
{analysis_name}_6_rgSet_vs_Noob_density_plots.pdf (genomewide density from RAW methylation values compared to after ssNoob is applied. Estimated from a random (set.seed) selection of 50000 CpGs.)

{analysis_name}_8_qc_plot.pdf (minfi qc plot - median intensity from methylated vs unmethylated intensities. A mean of these being less than 10.5 indicates bad sample/low signal)

{analysis_name}_11_sex_plot.pdf  - The median intensity signal from chromosome Y probes against X probes, indicating Male or Female
{analysis_name}_14_BMIQ_density_plots.pdf (similar to {analysis_name}_6_rgSet_vs_Noob_density_plots.pdf, but with genomewide density after BMIQ is applied as well)

{analysis_name}_15_histogram_of_differences.pdf (Comparison of RAW vs ssNoob methylation values. The difference of the two matrices are calculated, then the rowMeans (plot 1) and rowSds (plot 2) are calculated and plotted as histograms. The same is done between ssNoob methylation values and BMIQ methylation values. These plots are for documenting the difference after each major normalization method is applied.)

## tmp_results
tmp_results/:
(Moved to https://github.com/folkehelseinstituttet/mobagen/wiki/Methylation)
in this folder, all intermediate results are stored. These are not important for the end results, but can be inspected if the pipeline crashes at some point.

## results
results/:
{analysis_name}_2_SNP_Betas.rds (Matrix of beta values for the 65 (450k) or 60 (EPIC) SNPs on the array. Rows - rs id, column - samples)

{analysis_name}_3_control_probe_PCA.rds (The PCA object of PCA analysis from control probes)
{analysis_name}_9_estimated_cell_proportions.csv (Matrix with samples as rows and cell type proportions estimated as columns)
{analysis_name}_2_probes_removed.csv (2 column matrix of CpG-ids removed and why: "cross_hybridizing", "polymorphic" (SNP influenced))

{analysis_name}_11_added_sex_prediction_to_pheno.csv (Sample sheet data with additional columns for predicted sex, median intensity across Y chromosome probes, and median intenisty across X chromosome probes)

{analysis_name}_4_bisulphite_conversions.csv (The calculated bisulphite conversion rate for each sample)

{analysis_name}_4_num_probes_with_proportion_failed_samples_p0.01.csv (For each sample, the fraction of failed probes and the total of failed probes, compared to the detection p value cut off (0.01))
{analysis_name}_4_probes_failed_per_sample_p0.01.csv (For 5 different proportions (0.5, 0.1, 0.05, 0.01, 0.0033), how many CpGs have a bad detection p value for more than the give proportion of samples)

{analysis_name}_8_filtered_samples.csv (2 column matrix of sample ids removed and why: "NA_control_probes" (too many NA/INF values), "bisulphite" (low bisulphite conversion), "low_median_meth_or_unmeth_channel" (less than 10.5 mean median meth and unmeth log2-intensity), "manual_removed_bad_density" (removed due to having outlying values in either the control_probe_boxplot or genome-wide density), "pvalue" (removed samples with less than 90% detected probes based on detection pvalue cutoff))

# Quality Control Pipeline documentation/descrition

- load iDat files to RGset (all information in minfi object format)

- filter probes from RGset
	- poor quality probes (detection p-value larger than 0.01 for more than 5 % of the samples)
	- cross-hybridizing probes (probes shown to have more than 1 possible binding site in the genome)
	- non-CpG methylation (CpH, (H = A, C or T))
	- probes influenced by SNPs with higher than 0.05 MAF. SNP placed in interrogated CpG or at the neighboring positions

- create boxplot and PCA plot of control probes:
	- PCA pot: PC2 against PC1 from PCA of approximately 1200 Red and Green channel control probes
		- meant to use as documentation and to get intuition of possible amount of outliers
	- boxplot of control probe values for each sample
		- if a sample has outlying pattern, such as very narrow range of values -> remove sample

- filter samples with:
	- more than 25% missing/infinite values for control probes
	- more than 10% probes with high detection p value (> 0.01)
	- less than 90% bisulphite conversion rate (calculated with bscon() function, which uses bisulphite conversion control probes to estimate conversion rate)
	
- preprocess with ss-Noob:
	- single-sample-(normal-exponential out-of-band) method used for background subtraction and dye-bias correction. It uses out-of-band probes to estimate a normal-exponential background signal for each color channel, and normalizes the signal from regular probes based on these estimated distributions to obtain a purer signal where the background is removed. This algorithm only uses information within each sample to correct each sample. It has been shown to be one of the better background subtraction and dye-bias correction methods: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5408810/
	- from the authors of ssNoob: "We note that on the Beta value scale, there is no difference between values returned by Noob or ssNoob. Differences are confined to the methylated and unmethylated signals."

- plotting of genomewide densities of RAW methylation data and after ssNoob is applied
	- if a sample has clear outlying distribution from the others, remove sample

- removing samples based on outlying pattern

- check median methylated intensity against median unmethylated intensity, and remove samples with (median_M + median_U)/2 < 10.5. Indicates low signal from either one or both channels.

- cell type proportion estimation. Uses either Blood reference set or CordBlood reference set, dependening on input to config.yaml.
Blood reference set: DOI: 10.1371/journal.pone.0041361
CordBlood reference set: doi: 10.1080/15592294.2016.1161875

- sex prediction and removement of sex discordant samples between given sex and predicted sex
	- measures the median intensity values from combined methylated and unmethylated signal from Y chromsome probes and the same for X chromosome probes. If the differences between these medians is more than 2 (less than -2) a female is predicted.
	
- seperate out sex chromosome measurements from autosomal measurements

- do BMIQ: probe type normalisation. The Beadchips have 2 different probe designs, I and II. These have a clear difference in dynamic range for the obtained methylation measurements. BMIQ projects the global distribution of probe type II measurements onto the distribution of probe type I, meaning that the methylation values from probe type II measuements are slightly shifted. This is especially an usual step when analysing regions of data, i.e. looking for differentially methylation regions. It is also believed that the dynamic range of probe type I is more correct than probe type II. BMIQ also calculates only within each sample.

- plotting og genomewide distribution of methylation values after BMIQ is applied.

- plot histogram of differences between raw methylation values and after ssNoob is applied, and between methylation values after ssNoob is applied and BMIQ is applied. These plots are for documenting the difference after each major normalization method is applied.

- clean up and delete unnecassry data stored at intermediate steps


More information:

detection p value:
a detection p value greater than 0.01 should not be trusted. It is calculated by comparing the total dna signal (methylated + unmethylated) for each position and sample to the background signal level. The background is estimated using the negative control probes, assuming a normal distribution. Calculations are performed on the original (non-log) scale.

negative control probes: specifically design not to match the human genome, reads/signal from these probes in red and green channel to estimate the background signal distribution.

out of band probes: for a given type I bead, the intensity from the unused color channel has been proposed as a means of estimating background signal, and termed the out-of-band intensity. 



