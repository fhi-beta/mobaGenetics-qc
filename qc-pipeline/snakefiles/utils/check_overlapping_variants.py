import pysam
import pandas as pd
from pathlib import Path
from collections import defaultdict

def process_vcf_files(batch_genotyping_chip, chromosomes, tmpMod6_path, output_file):
    with open(output_file, 'w') as f:
        for chromosome in chromosomes:
            print(f"Processing chromosome {chromosome}...")
            variants_by_batch = defaultdict(set)
            for batch in batch_genotyping_chip.keys():
                print(f"Reading batch {batch}...")
                vcf_path = tmpMod6_path / batch / f"mod6_conform.chr{chromosome}.conformed.bcf"
                with pysam.VariantFile(str(vcf_path), 'r') as vcf:
                    for record in vcf:
                        variants_by_batch[batch].add((record.chrom, record.pos, record.ref, record.alts))

            df = pd.DataFrame(index=batch_genotyping_chip.keys(), columns=batch_genotyping_chip.keys(), data=0)

            for batch1 in batch_genotyping_chip.keys():
                for batch2 in batch_genotyping_chip.keys():
                    if batch1 == batch2:
                        df.loc[batch1, batch2] = len(variants_by_batch[batch1])
                    else:
                        overlap = variants_by_batch[batch1].intersection(variants_by_batch[batch2])
                        df.loc[batch1, batch2] = len(overlap)

            f.write(f"# Overlapping Variants for Chromosome {chromosome}\n\n")
            f.write(df.to_markdown())
            f.write("\n\n")


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
    "snp014": "Illumina Global Screening Array MD v1.0",
    "snp015a": "DeCodeGenetics 1_v2",
    "snp015b": "DeCodeGenetics v3 1_v2",
    "snp016a": "DeCodeGenetics v3 1_v2",
    "snp016b": "DeCodeGenetics v3 1_v2",
    "snp017a": "DeCodeGenetics v3 1_v2",
    "snp017b": "DeCodeGenetics v3 1_v2",
    "snp017c": "DeCodeGenetics v3 1_v2",
    "snp017d": "DeCodeGenetics v3 1_v2",
    "snp017e": "DeCodeGenetics v3 1_v2",
    "snp017f": "DeCodeGenetics v3 1_v2",
    "snp018a": "DeCodeGenetics v3 1_v2",
    "snp018b": "DeCodeGenetics v3 1_v2",
    "snp018de": "DeCodeGenetics v3 1_v2"
}
chromosomes = range(1, 23)
tmpMod6_path = Path("/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod6-imputation/")
output_file = "/mnt/work/oystein/tmp/variant_overlap_tables.md"

process_vcf_files(batch_genotyping_chip, chromosomes, tmpMod6_path, output_file)