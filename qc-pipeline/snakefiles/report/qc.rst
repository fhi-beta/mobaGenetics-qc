(A demo to short what could have been included in this section)

QC report for {{ snakemake.config["dataSet"] }}
Chip used {{ snakemake.config["chip"] }}

For a general description of the QC pipeline, see https://github.com/folkehelseinstituttet/mobaGenetics-qc/wiki (this link does not work for all yet as the wiki is not public)

This report has two kind of Results_ (click the header to expand the result section)

The ones that describes the output main modules:

* `- Module 1 Final report`_

The ones that decribes the intermediate results within a module 

* `Module 1 Data preparation`_
* `Module 2 Core samples and infere pedigree`_

And here is the link to a plot `maf_removal_markers.png`_


Temporary example that show that we could show parameters here.
Cut off parameters (not updated yet):

====================  ====================
Parameter                   Value 
====================  ====================
Cluster Separation       < {{ snakemake.config["cluster_sep_thr"] }}
10% GC Score             < {{ snakemake.config["10%_GC_score_thr"] }}
AA theda dev             > {{ snakemake.config["aa_theta_dev_thr"] }}
====================  ====================


