#### Table of contents <!-- :TOC: -->
- [Introduction](#introduction)
  - [Planned QC stages](#planned-qc-stages)
- [Setup](#setup)
  - [Software install](#software-install)
  - [Configure](#configure)
  - [Run](#run)
  - [Results](#results)
- [Todo](#todo)

# Introduction
( when Gutorm edits this in emacs it creates the nifty table of content automagically.)

## Planned QC stages
The overall structure of the QC will probably have the following stages

1. Prepare a curated public raw-data set with bedset having correct
   .fam file with respect to parents and sex. (This is not completely
   done yet). There will be one such bedset for each data-set. 
1. Create a QC-pipline that can be run on all sets individually and
   that ultimately creates pre-imputation bedset. This work is in
   progress as of June 2022 and most of the pipeline has been written. 
1. Merge all sets that have been genotyped with the same chip, but
   also remove duplicate individuals within the merege sets. This will
   be called the chipMerged set and will be imputeded as well.
1. Merge all chipMerged imputation results to a final imputed
   set. Duplicate individuals will again be removed.


# Setup

## Software install
Clone the git repository. 

Then install (if not installed)
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
An exported environment (suitable for running the pipeling) is
available on `mobaGenetics-qc/aux/environment/snpQc-env.txt` . The
environemnt is suitable for a linux-64 platform and can be used to
create an environment (provided `conda` is installed) by

`conda create --name snpQc --file snpQc-env.txt`

Among many things it will install R, python and snakemake for you. When done, 
acitavte it:

`conda active snpQc`

and you should be good to go. 

### Nitty gritty details on the conda install

All packages are grabbed form conda-forge, it might be we change this
the above to using mamba and mamba-forge.

The `snpQc-env.txt' was created by the
command

`conda list --explicit > snpQc-env.txt`

in an (ubuntu) environemnt (132Gb) that successfully runs the QC. In
theory you still could get problems, let us know if you identify some - like  minimum
memory requirements.

A yaml version might be produced later, it would be more robust so it
would work on say Windoies.  For now we will not maintain two
configuration files. If you want a yaml-version you can always do 

`conda env export > snpQc-env.yml`

The yaml-file should the be used to recreate the environment with 

`conda env create --file snpQc-env.yml`

## Configure 
You need to edit the `snakefiles/config.yaml` file until it works (many paths need editing)

## Run

`snakemake -c8 rulename`

As of 29.6.2022 you have to run 

`snakemake -c8 rayner_report` and then `snakemake -c8` for the
remaining. When the default rule is fixed, the first command will not
be necessary.

should work but does not yet work for all rules. (Due to certain steps
not being finished and the default rule being poorly setup)

## Results
Check the directory `output_base` from `config.yaml`

The `results` subdirectory contais the results of (almost?) every
rule. The special config file `rules.yaml` is used to document the
rules and is parsed to create a propper flow.

For example, `plink` commands in the pipeline will typically take an input bedset,
create an output bedset and some python code will summarize the difference
between the two and put it in the `results` directory (see below). 

The format should be shuch that flag-files (or interactive queries) can be made here - a not finnished script 

A html-report can be generated by 

`snakemake --report reportname.html rule`

after `rule` being previously run. It will create a html-report
described below (default name being report, but you typically want the
set-name here, sucha as 'snp014.htm')

The subdirectory `report` (same level as Snakefile) contains files to
that are inserted in the big html-file - picked up from `results` -
see the section below for more details.

### Html-report

The QC report (see how to make above) is produced by [standard
Snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html)
. It is itself a product under development, so functionality might
change or might not be optimally used.

It is not obvious that this will be the final report, but we are
piloting it. Also keep in mind that since the code was written in
2019/2020 - new functionality might has been added in snakemake that
we are not using.

There is a master file for the QC-report residing in the
snakefiles/report/qc.rst . This files contains the first part of the
reports explaining what the report is.

What goes into the report is configured through snakemake `output:`
elements withing the `report()`command - typically one or more for
each rules.

Refer to the snakemake documentation for syntax etc

In order for the report to have entries sorted, the results section
use a description field that can be sorted. The first line of the
itself is configured in the field `Rule order` in the `rules.yaml`
found with the snakefiles.  The snakemake files (modx*) will, for each rules
that has a report, create a corresponding .rst file. This contains a
sortable header (see over) as well as a recap of the rule results.

`rules.yaml`also contain easilly editable documentation for the rules
that will be picked up and inserted in the html-file.

Note that this is a way of cramming results into the html-report,
making it easy to quickly browse what the rule did. It is also a
stretch of the functionality proposed by snakemake ... More details
can be found in yaml-files (the File column of the report).

In addition to the .yaml file, the pipeline often produces as
.yaml.details file. The latter contains detailed information like what
samples/markers that were removed by the rule. The .details file is
not included, as all results files are embedded in the html-files,
producing a big file (typically 5-10 Mb, depending on the number of
plots).


### Intermediate results
Each module will have a separate directory where all the intermediate
results are stored. These are typically `plink` bedsets for every plink
command issued, including log files for debugging. 

On the directory `aux/bin`you will find an early pilot
(`createFlagList.py`) that take as input results (typically
`*.details` files generated from the pipeline) and creating flaglists
(where did samples fail/disapear) based on a template file. An example
of such a template file is found on `flagTemplate.csv`.

The flagfiles is just one way of visualizing what has happened during
the pipeline run - this was done for MobaGenetics 1.0. It might be
that the new pipeline can offer the same information in other forms -
an example would be a script that takes one or more samples as input
and show everything that happened to that sample during the pipeline.


# Todo
Make subsections here to visualize bulks of work and priotizing. 
- Change the config-file:
  - one global (paths, tresholds etc)
  - one local (name of the chip, overrun of the global values)
- Split out the parts of the pipeline that 
  - Rename ids (should come from the pub directory)
  - Set parents (should come from the pub directory as well)
  - Set sex (should be right in the pub directory)
- Checks & and repairs pedigree (should be right in the pub directory)
- Clean up directories in the source, files are now all over the place (aux/bin/scripts/docs)
- The rule `common_markers_moba_ref`: takes forever. The reduction of
  the 1000 genomes set. It could be that this one should be isolated,
  right now there is a terrible clumsy 'uncomment this part unless you
  want this to take forever'. See relevant code in mod3.

- Clean up the default rule so they actually run to the end of the pipeline. 

- Gutorm also had problems in january 2020 in mod3 after
  `common_markers_moba_ref` where data were projected on 1000g.There
  was problem with the plotting itself (could have been memory issues
  on Harvest at the time) but generally speaking, something is rotten
  using the 1000g reference. For the actuall rotterdam2 pipeline
  hapmap was used.

- Output
  - eigenvalues/vectors/pve end up on the install-directory, not the results directory
