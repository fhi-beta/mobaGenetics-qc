(This html file is autogenrated by snakemake)

QC report for {{ snakemake.config["dataSet"] }}.
Chip used was {{ snakemake.config["chip"] }}

For a general description of the QC pipeline, see
https://github.com/folkehelseinstituttet/mobaDocsBackoffice/wiki/pipeline
(this link does not work for all yet as the wiki is not public)

This report has two kind of Results (click the header of the Result
section in the sidebar to expand the result section)

The ones that describes the output main modules:

* `- Module 1 Data conversion recap`_
* `- Module 2 Pedigree fix recap`_
* `- Module 3 Good markers recap`_
* `- Module 4 Core individuals recap`_
* Module 5: Not ready yet


The ones that decribes the intermediate results within a module 

* `Module 1 Data preparation`_
* `Module 2 Core samples and infere pedigree`_
* `Module 3 Good markers`_
* `Module 4 Unrelated samples`_
* `Module 5 Phasing and imputation preparation`_


The result column *Description* can be used to sort the rules in the
order they appear in the pipleine. The *File* column can be used to
see more details from the rules.

Sandbox: Some examples of other types of information:
    
Links: A link to a plot: `maf_removal_markers.png`_ 
Table: (note that these values are shown in the result section anyway, this is just a demo of functionality that could have been included in the report).

====================  ====================
Parameter                   Value 
====================  ====================
Cluster Separation       < {{ snakemake.config["cluster_sep_thr"] }}
10% GC Score             < {{ snakemake.config["10%_GC_score_thr"] }}
AA theda dev             > {{ snakemake.config["aa_theta_dev_thr"] }}
====================  ====================


