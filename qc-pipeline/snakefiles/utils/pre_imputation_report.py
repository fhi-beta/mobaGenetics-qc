import pandas as pd
def write_report(output_filename, batch, bedset):
    md_file = open(output_filename, "a")
    md_file.write(f"#Pre-imputation report for batch {batch}")
    fam_df = pd.read_csv(bedset[2], delim_whitespace=True, header=None,
    names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])
    included_number_of_samples = fam_df.shape[0]
    md_file.write(f"\n{included_number_of_samples} samples included")
    md_file.close()