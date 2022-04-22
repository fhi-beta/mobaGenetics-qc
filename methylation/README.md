#### Table of contents <!-- :TOC: -->
- [Introduction](#introduction)
- [Setup](#setup)
  - [TSD](#tsd)
- [Configuration](#configuration)
  - [Global config file: globalConfig.yaml](#global-config-file-globalconfigyaml)
  - [Local config-files](#local-config-files)
- [Pre-quality control:](#pre-quality-control)
- [Running the pipeline](#running-the-pipeline)
  - [Canned wrapper script](#canned-wrapper-script)
  - [Removing bad samples](#removing-bad-samples)
  - [Finalizing the run](#finalizing-the-run)
- [Results](#results)
  - [Final](#final)
  - [Logs](#logs)
  - [Plots](#plots)
  - [tmp_results](#tmp_results)
  - [results](#results-1)
- [Quality Control Pipeline documentation/description](#quality-control-pipeline-documentationdescription)

# Introduction
Based on prior QC work by Christian M. Page and Haakon E. Nustad and
the bhklab/IlluminaEPICmethylation pipeline on github by Christopher Eeles and Joanna Pryzbyl.
Developed by Hanna Zdanowicz, Haakon E. Nustad and Gutorm  Høgåsen.

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
- output_path is where all the results logs and whatnot will end
  up. The pipeline might create terrabytes of data depending on your
  input
- max_sample_size_* Unless you have loads of memory, you want to
  reduce the size of samples being processed in one go. See the
  documentation and comments in the yaml-file for more details. If
  your samples size within a set is larger than this amount, the
  pipeline will return an error message the first time it runs, but
  will create temporary config files where the data is split into
  suitable batches.  This splitting process is entirely handled by the
  [Canned wrapper script](#canned-wrapper-script) described below
- renv If you cannot use conda as described below, an alternative
  using pip & renv is possible.
  

## Local config-files 
Each dataset has one or more config-file (called a local config file)
that can be used to override the global one. It must be on the form
setname.yaml. Sometime the setname will include child or parent. 

If you use the standard/original 'MoBa' path-structure you don't need
to change anything.

Important variables:

- set\_name: The name of the set, will be used among other things to
  name the directories under output_path found in the global config
  file described above.
- data\_type: Chip name
- cell\_type: We assume all samples in the set have the same type.
- local\_path: Where to find idat files
- discarded\_samples: List of samples to be ignored, typically results
  of manual inspection. See [removing bad samples](#Removing bad
  samples) below.

Go through the local config.yaml file (usual names met003.yaml,
met004.yaml, etc) in input folder and change parameters to your
specific file paths and names, especially data directory (idat and
sample sheet file path), data type (450k or EPIC). 

Additionally, add a file in input folder:
manual\_remove\_samples\_bad_density.txt. This file is meant to contain
sample ids of samples you want to remove based on having outlying
genomewide DNA methylation patterns. This should be filled with sample
ids based on output from the first and second part of the
pipeline. Good to go!


# Pre-quality control:
Important pre-step: the sample sheet needs the beadchip array column
to have column name: Sentrix\_ID, and the position on the array to have
column name: Sentrix\_Position. It also needs a column sampleType (see
below) This is important for minfi to automatically recognize correct
columns to identify the .idat files.

All samples you wish to run QC for must be listed in your
sample\_sheet.csv, which must be located in the folder you specify in
config.yaml. It should be placed together with the .idat files. 
This .csv file is recognized automatically with the minfi function applied
to read in the sample_sheet. 

The sample sheet must contain a column SampleType in order to detect
if you have a mix of children and older individuals by the script
0\_checkSampleSize.R (snakemake rule sample\_size\_check). If so, the script
will create an error with an error message, encouraging the user to
split their samples into two .csv files, and running these
seperately. This is because the quality of the QC will increase if
such biologically related differences are not present.

# Running the pipeline
The first time the pipeline is run, remove bad samples manually as
described below by using the (snakemake) rules first\_part and second_part.

In order to run Snakemake, most users will want to use the [Canned wrapper
script](#canned-wrapper-script) described below. 

The sections further below describe the different partst of the
pipleline and how to run them with snakemake commands.


## Canned wrapper script 
The wrapper script runOneorMore.sh will make the process of running
the pipeline easier, especially if your resources are limited. 

At the base you will be running commands like 

$ snakemake --core 1
--configfile input/globalConfig.yaml met001.yaml --snakefile
Snakefile third\_part

where the last parameter is the snakemake rule (workflow part) you
want to run. The canned script echoes what it does, so you don't need
to read up on [snakemake](https://snakemake.readthedocs.io/en/stable/)
(for now).

Do not create the final results (rule third\_part) before bad samples
have been removed, see the section below that describes how to create
setname\_manual\_remove\_samples\_bad\_density.txt files.


You will need to edit the script runOneorMore.sh to set some crucial
variables:

- output\_path Where you want the results
- data\_path Where idat-files are found. For MoMics users this will be
  .../Methylation/Datasets
- global\_config Where to find the globalConfig.yaml described in a
  separate section. It sets common parameters for all datasets you
  will run. For MoMics users this is set right as long as you have set
  data\_path correctly (see above).
- rule What part of the pipeline you want to run This should in most
  cases be third_part. The exception is the first time you run the
  pipeline as you need to generate lists over bad samples - for MoMics
  users these are already found and placed under the sets QC/input
  sub-directory.
- cores Max number of cores to use

If you wish to run one local config file, met001.yaml for example, you
do this with the following command:

bash runOneorMore.sh path/met001.yaml 

For MoMics users the path will be .../Methylation/Datasets/met001/QC/input/

You can give several input config files if desired (the syntax assumes
bash and will run met001, met004 and met007):

bash runOneorMore.sh Datasets/met\*/QC/input/met00[147]\*.yaml 

If the amount of samples in a set are more than the max\_sample\_size
specified in globalConfig.yaml, the data will be split automatically
into met001-1, met001-2, met001-3, etc. and run for each if these
batches. 

### A warning on failing rule sample\_size\_check
When snakemake is run (usually by runOneorMore.sh) now rules should
fail _except_ sample\_size\_check.

This rule detects that there are two many samples in the set to handle
for your limited memory. It will split the set up using a directory
tmp_config. The results will be a lot of subsets that Snakemake will
be run on.

## Removing bad samples
The pipeline is organized in 3 parts, where one should inspect the
plots produced located in plots/ after the first and second part. This
section describes the two first parts, leading to boxplots suitable
for removing bad samples. 

Based on the plots produced in first and second part, you should add
sample IDs to the file input/(x)\_manual\_remove\_samples\_bad\_density.txt
if you think they should be removed. Remove samples showing especially
outlying patterns compared to the rest, based on the boxplot of
control\_probes and genome wide density plot.


## Finalizing the run

After one has run the third\_part, the following command will clean up
the output files, removing unnecessary files:

$ snakemake --core 1 --configfile input/globalConfig.yaml input/met001.yaml --snakefile Snakefile clean\_up


# Results
All results are stored within the output\_path folder (see above), in
subdirectories named after the dataset.  (config file, minus the .yaml
part) This is set automatically.

You will notice that some sets will be split (due to memory limitation
that you can tune in the in global configfile, see above.

## Final 
The end matrices with data:

This section has been moved completely to https://github.com/folkehelseinstituttet/mobagen/wiki/Methylation#QC

## Logs
Additional files with valuable information:
logs/:
- The Rout files are log-files containing information about the R
  session and the functions applied to the data. These will be found
  on the logs subdirectory where the Results are.  They contain
  information of what is done to the data during the different
  scripts. Additionally, if an error occurs in a script, the error
  message is given in the associated log file. These files are
  primarily meant for debugging the pipeline, and can be inspected
  to understand the order of when the different functions and QC steps
  have been applied.

## Plots
plots/:
(Moved to https://github.com/folkehelseinstituttet/mobagen/wiki/Methylation)

## tmp_results
tmp_results/:

in this folder, all intermediate results are stored. These are not important for the end results, but can be inspected if the pipeline crashes at some point.

## results
Section moved to
https://github.com/folkehelseinstituttet/mobagen/wiki/Methylation 

# Quality Control Pipeline documentation/description
Section moved to
https://github.com/folkehelseinstituttet/mobagen/wiki/Methylation 
