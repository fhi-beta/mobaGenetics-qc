QC report for {{ snakemake.config["dataSet"] }}

Cut off paramenters:


====================  ====================
Parameter                   Value 
====================  ====================
Cluster Separation       < {{ snakemake.config["cluster_sep_thr"] }}
10% GC Score             < {{ snakemake.config["10%GC_score_thr"] }}
AA theda dev             > {{ snakemake.config["aa_theta_dev_thr"] }}
====================  ====================


