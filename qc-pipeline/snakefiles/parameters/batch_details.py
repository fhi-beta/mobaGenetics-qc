#!/usr/bin/python3

batch_genotyping_chip = {
    "snp001": "Illumina Human Core Exome 12 v.1.1",
    "snp002": "Illumina Human Core Exome 12 v.1.1",
    "snp003": "Illumina Human Core Exome 24 v.1.0",
    "snp012": "Illumina Global Screening Array MD v1.0",
    "snp014": "Illumina Global Screening Array MD v1.0"
}
batch_mobagenetics_10_name = {
    "snp001": "Harvest m12b",
    "snp002": "Harvest m12a",
    "snp003": "Harvest m24",
    "snp012": "Rotterdam1",
    "snp014": "Rotterdam2"
}
batch_snp_table = {
    "snp012": "/mnt/archive2/MomicsSource/snpArray/snp012/ROTTERDAM1/rotterdam1-aux/qc/snp_table_rotterdam1.txt",
    "snp014": "/mnt/archive2/MomicsSource/snpArray/snp014/ROTTERDAM2/delivery-fhi/data/raw-data/snptable/snp_table_rotterdam2.txt"
}

# Returns the chip of a given batch
def getChip(batch):
    return batch_genotyping_chip[batch]

# Returns the snp table of a given batch
def getSnpTable(batch):
    return batch_snp_table[batch]


