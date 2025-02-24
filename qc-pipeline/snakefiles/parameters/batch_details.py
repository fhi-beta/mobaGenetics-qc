#!/usr/bin/python3

batch_genotyping_chip = {
    "snp001": "Illumina Human Core Exome 12 v.1.1",
    "snp002": "Illumina Human Core Exome 12 v.1.1",
    "snp003": "Illumina Human Core Exome 24 v.1.0",
    "snp007": "Illumina Human Omni Express 24 v.1.0",
    "snp008": "Illumina Human Omni Express 24 v.1.0",
    "snp009": "Illumina Infinium Omni Express 24 v.1.2",
    "snp010": "Illumina Global Screening Array MD v1.0", # Needs to be checked
    "snp011": "Illumina Infinium Omni Express 24 v.1.2",
    "snp012": "Illumina Global Screening Array MD v1.0",
    "snp013": "DeCodeGenetics V1_v2", # Needs to be checked
    "snp014": "Illumina Global Screening Array MD v1.0",
    "snp015a": "Illumina Global Screening Array MD v1.0",
    "snp015b": "Illumina Global Screening Array MD v3.0",
    "snp016a": "Illumina Global Screening Array MD v3.0",
    "snp016b": "Illumina Global Screening Array MD v3.0",
    "snp017a": "DeCodeGenetics v3 1_v2",
    "snp017b": "DeCodeGenetics v3 1_v2",
    "snp017c": "DeCodeGenetics v3 1_v2",
    "snp017d": "DeCodeGenetics v3 1_v2",
    "snp017e": "DeCodeGenetics v3 1_v2",
    "snp017f": "DeCodeGenetics v3 1_v2",
    "snp018a": "DeCodeGenetics v3 1_v2",
    "snp018b": "DeCodeGenetics v3 1_v2",
    "snp018c": "DeCodeGenetics v3 1_v2",
    "snp018de": "DeCodeGenetics v3 1_v2",
    "snp019": "GSAMD-24v2-0_20024620_B1"
}
batch_name = {
    "snp001": "Harvest m12b",
    "snp002": "Harvest m12a",
    "snp003": "Harvest m24",
    "snp007": "Norment Jan15",
    "snp008": "Norment Jun15",
    "snp009": "Norment May16",
    "snp010": "Norment Feb18",
    "snp011": "Ted1",
    "snp012": "Rotterdam1",
    "snp013": "Ted2",
    "snp014": "Rotterdam2",
    "snp015a": "Norment Feb20 1",
    "snp015b": "Norment Feb20 3",
    "snp016a": "Norment Aug20 996",
    "snp016b": "Norment Aug20 1029",
    "snp017a": "DeCode 1",
    "snp017b": "DeCode 1",
    "snp017c": "DeCode 1",
    "snp017d": "DeCode 1",
    "snp017e": "DeCode 1",
    "snp017f": "DeCode 1",
    "snp018a": "DeCode 1",
    "snp018b": "DeCode 1",
    "snp018c": "DeCode 1",
    "snp018de": "DeCode 1",
    "snp019": "LifeBrain"
}

batch_snp_table = {
    "snp011": "/mnt/archive/Momics/MomicsPub/snpArray/Datasets/snp011/raw-data/aux/snp-table.txt",
    "snp012": "/mnt/archive/Momics/MomicsPub/snpArray/Datasets/snp012/raw-data/aux/snp_table.txt",
    "snp014": "/mnt/archive/Momics/MomicsPub/snpArray/Datasets/snp014/raw-data/aux/snp_table.txt"
}

batch_cluster_sep_col_name = {
    "snp011": "InfiniumOmniExpress-24v1-2_A1.bpm.Cluster Sep",
    "snp012": "GSAMD-24v1-0_20011747_A4.bpm.Cluster Sep",
    "snp014": "GSAMD-24v1-0_20011747_A4.bpm.Cluster Sep"
}

batch_aa_theta_dev_col_name = {
    "snp011": "InfiniumOmniExpress-24v1-2_A1.bpm.AA T Dev",
    "snp012": "GSAMD-24v1-0_20011747_A4.bpm.AA T Dev",
    "snp014": "GSAMD-24v1-0_20011747_A4.bpm.AA T Dev"
}

batch_founders_offspring = {
    "snp007": "founders",
    "snp008": "founders",
    "snp018b": "founders",
    "snp018de": "founders"
}

batch_ab_allele_mapping_files = {
    "snp001": "ab_conversion/HumanCoreExome-12v1-1_B.update_alleles.txt.gz",
    "snp002": "ab_conversion/HumanCoreExome-12v1-1_B.update_alleles.txt.gz",
    "snp003": "ab_conversion/HumanCoreExome-24v1-0_A.update_alleles.txt.gz"
}

batch_strand_mapping_files = {
    "snp001": "ab_conversion/HumanCoreExome-12v1-1_B-b37.Ilmn_chr.strand.gz",
    "snp002": "ab_conversion/HumanCoreExome-12v1-1_B-b37.Ilmn_chr.strand.gz",
    "snp003": "ab_conversion/HumanCoreExome-24v1-0_A-b37.Ilmn_chr.strand.gz"
}

batch_exclusion_files = {
    "snp001": ["exclusion_harvest/CHARGE_ExomeChip_v1.0_Excluded_Variants.txt", "exclusion_harvest/pchip_blackList_dec2015_stripped.txt"],
    "snp002": ["exclusion_harvest/CHARGE_ExomeChip_v1.0_Excluded_Variants.txt", "exclusion_harvest/pchip_blackList_dec2015_stripped.txt"],
    "snp003": ["exclusion_harvest/CHARGE_ExomeChip_v1.0_Excluded_Variants.txt", "exclusion_harvest/pchip_blackList_dec2015_stripped.txt"]
}

initial_batch_merging = {
    "snp018de": ["snp018d", "snp018e"]
}

# Returns the chip of a given batch
def getInitialMerging(batch):
    if batch in initial_batch_merging:
        return initial_batch_merging[batch]
    else:
        return [batch]

# Returns the file needed for AB allele mapping, "None" if no mapping is needed.
def getAbAlleleMappingFile(batch):
    return batch_ab_allele_mapping_files.get(batch, "None")

# Returns the file needed for strand mapping, "None" if no mapping is needed.
def getStrandMappingFile(batch):
    return batch_strand_mapping_files.get(batch, "None")

# Returns a list of variant exclusion files
def getBlackListFiles(batch):
    return batch_exclusion_files.get(batch, [])

# Indicates whether only founders or both founders and offspring need to be QCed
def getFoundersOffspring(batch):
    if batch in batch_founders_offspring:
        return [batch_founders_offspring[batch]]
    else:
        return ["founders", "offspring"]

# Returns the chip of a given batch
def getChip(batch):
    return batch_genotyping_chip[batch]

# Returns the snp table of a given batch
def getSnpTable(batch):
    return batch_snp_table.get(batch, "utils/dummy")

# Returns the column for cluster separation in a given snp table
def getClusterSeparationColumn(batch):
    return batch_cluster_sep_col_name.get(batch, "None")

# Returns the column for cluster separation in a given snp table
def getAaThetaDevColumn(batch):
    return batch_aa_theta_dev_col_name.get(batch, "None")


