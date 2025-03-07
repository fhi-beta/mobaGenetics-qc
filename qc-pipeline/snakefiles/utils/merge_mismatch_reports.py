
import pandas as pd

mismatch_information_list = []
mismatch_relationship_list = []

batches = ['snp001', 'snp002', 'snp003', 'snp007', 'snp008', 'snp009', 'snp010', 'snp011', 'snp012', 'snp014', 'snp015a', 'snp015b', 'snp016a', 'snp016b', 'snp017a', 'snp017b', 'snp017c', 'snp017d', 'snp017e', 'snp017f', 'snp018a', 'snp018b', 'snp018c', 'snp018de']

for batch in batches:
    mismatch_information_path = '/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod2-genetic-relationship' /batch /"mismatch_information.gz"
    mismatch_information_table = pd.read_csv(mismatch_information_path, delimiter = '\t')
    mismatch_information_table["batch"] = batch
    mismatch_information_list.append(mismatch_information_table)

    mismatch_relationship_path = '/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod2-genetic-relationship' /batch /"mismatch_relationship.gz"
    mismatch_relationship_table = pd.read_csv(mismatch_relationship_path, delimiter = '\t')
    mismatch_relationship_table["batch"] = batch
    mismatch_relationship_list.append(mismatch_relationship_table)



mismatch_information_table = pd.concat(mismatch_information_list, ignore_index = True)
mismatch_relationship_table = pd.concat(mismatch_relationship_list, ignore_index = True)

mismatch_information_table.to_csv('/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod2-genetic-relationship/mismatch_information.gz', sep = '\t', index = False, compression = 'gzip')
mismatch_relationship_table.to_csv('/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod2-genetic-relationship/mismatch_relationship.gz', sep = '\t', index = False, compression = 'gzip')

