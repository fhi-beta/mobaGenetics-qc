#!/usr/bin/python3

batch_genotyping_chip = {
    "snp001": "Illumina Human Core Exome 12 v.1.1",
    "snp002": "Illumina Human Core Exome 12 v.1.1",
    "snp003": "Illumina Human Core Exome 24 v.1.0",
    "snp007": "Illumina Human Omni Express 24 v.1.0",
    "snp008": "Illumina Human Omni Express 24 v.1.0",
    "snp009": "Illumina Infinium Omni Express 24 v.1.2",
    "snp011": "Illumina Infinium Omni Express 24 v.1.2",
    "snp012": "Illumina Global Screening Array MD v1.0",
    "snp014": "Illumina Global Screening Array MD v1.0",
    "snp015a": "Illumina Global Screening Array MD v1.0",
    "snp015b": "Illumina Global Screening Array MD v3.0",
    "snp016a": "Illumina Global Screening Array MD v3.0",
    "snp016b": "Illumina Global Screening Array MD v3.0"
}
batch_name = {
    "snp001": "Harvest m12b",
    "snp002": "Harvest m12a",
    "snp003": "Harvest m24",
    "snp007": "Norment Jan15",
    "snp008": "Norment Jun15",
    "snp009": "Norment May16",
    "snp011": "Ted1",
    "snp012": "Rotterdam1",
    "snp014": "Rotterdam2",
    "snp015a": "Norment Feb20 1",
    "snp015b": "Norment Feb20 3",
    "snp016a": "Norment Aug20 996",
    "snp016b": "Norment Aug20 1029"
}

batch_snp_table = {
    "snp011": "/mnt/archive2/MomicsSource/snpArray/snp011/TED/delivery-fhi/data/raw-data/snp-table/snp-table.txt.gz",
    "snp012": "/mnt/archive2/MomicsSource/snpArray/snp012/ROTTERDAM1/rotterdam1-aux/qc/snp_table_rotterdam1.txt",
    "snp014": "/mnt/archive2/MomicsSource/snpArray/snp014/ROTTERDAM2/delivery-fhi/data/raw-data/snptable/snp_table_rotterdam2.txt"
}

batch_cluster_sep_col_name = {
    "snp011": "InfiniumOmniExpress-24v1-2_A1.bpm.Cluster Sep",
    "snp012": "GSAMD-24v1-0_20011747_A4.bpm.Cluster Sep",
    "snp014": "GSAMD-24v1-0_20011747_A4.bpm.Cluster Sep"
}

batch_aa_theta_dev_col_name = {
    "snp011": "GSAMD-24v1-0_20011747_A4.bpm.AA T Dev",
    "snp012": "GSAMD-24v1-0_20011747_A4.bpm.AA T Dev",
    "snp014": "InfiniumOmniExpress-24v1-2_A1.bpm.AA T Dev"
}

# Returns the chip of a given batch
def getChip(batch):
    return batch_genotyping_chip[batch]

# Returns the snp table of a given batch
def getSnpTable(batch):
    if batch in batch_snp_table:
        return batch_snp_table[batch]
    else:
        return "utils/mod1/dummy_snp_table"

# Returns the column for cluster separation in a given snp table
def getClusterSeparationColumn(batch):
    if batch in batch_cluster_sep_col_name:
        return batch_cluster_sep_col_name[batch]
    else:
        return "cluster_separation"

# Returns the column for cluster separation in a given snp table
def getAaThetaDevColumn(batch):
    if batch in batch_aa_theta_dev_col_name:
        return batch_aa_theta_dev_col_name[batch]
    else:
        return "theta_dev"


