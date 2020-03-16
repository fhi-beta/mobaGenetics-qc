(A demo to short what could have been included in this section)

QC report for {{ snakemake.config["dataSet"] }}
Chip used {{ snakemake.config["chip"] }}

For a general description of the QC pipeline, see https://github.com/folkehelseinstituttet/mobaGenetics-qc/wiki (this link does not work for all yet as the wiki is not public)

This report has two kind of Results_ (click the header to expand the result section)

The ones that describes the output main modules:

* `- Module 1 Data conversion recap`_
* `- Module 2 Pedigree fix recap`_
* `- Module 3 Good markers recap`_
* `- Module 4 Core individuals recap`_
* Not ready yet


The ones that decribes the intermediate results within a module 

* `Module 1 Data preparation`_
* `Module 2 Core samples and infere pedigree`_
* `Module 3 Good markers`_
* `Module 4 Unrelated samples`_
* `Module 5 Phasing and imputation preparation`_


 Some examples of other types of information:
  
  
Links: A link to a plot: `maf_removal_markers.png`_ 


Table: (note that these values are shown in the result section anyway) ...

====================  ====================
Parameter                   Value 
====================  ====================
Cluster Separation       < {{ snakemake.config["cluster_sep_thr"] }}
10% GC Score             < {{ snakemake.config["10%_GC_score_thr"] }}
AA theda dev             > {{ snakemake.config["aa_theta_dev_thr"] }}
====================  ====================


