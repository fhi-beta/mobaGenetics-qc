import pandas as pd
def write_report(output_filename, batch, bedset):
    md_file = open(output_filename, "a")
    md_file.write(f"# Pre-imputation report for batch {batch}")
    fam_df = pd.read_csv(bedset[2], delim_whitespace=True, header=None, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])
    md_file.write(f"\n## Samples overview")
    included_number_of_samples = fam_df.shape[0]
    md_file.write(f"\n{included_number_of_samples} samples")
    has_mother=fam_df[fam_df["MID"] != "0"]
    n_has_mother = has_mother.shape[0]
    has_mother_in_data = fam_df[(fam_df["MID"] != "0") & (fam_df["MID"].isin(fam_df["IID"]))]
    n_has_mother_in_data = has_mother_in_data.shape[0]
    n_missing_mothers=n_has_mother-n_has_mother_in_data
    has_father=fam_df[fam_df["FID"] != "0"]
    n_has_father = has_father.shape[0]
    has_father_in_data = fam_df[(fam_df["FID"] != "0") & (fam_df["FID"].isin(fam_df["IID"]))]
    n_has_father_in_data = has_father_in_data.shape[0]
    n_missing_fathers=n_has_father-n_has_father_in_data
    md_file.write(f"\n{n_has_mother} offspring with mother ID")
    md_file.write(f"\n{n_has_mother_in_data} offspring with mother ID in dataset")
    md_file.write(f"\n{n_missing_mothers} mothers missing from dataset")
    md_file.write(f"\n{n_has_father} offspring with father ID")
    md_file.write(f"\n{n_has_father_in_data} offspring with father ID in dataset")
    md_file.write(f"\n{n_missing_fathers} fathers missing from dataset")
    md_file.close()