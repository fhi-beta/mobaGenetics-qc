import pandas as pd
import matplotlib.pyplot as plt
import mobaQcTools as mqc
def write_report(output_filename, batch, bedset, imiss, lmiss):
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
    has_father=fam_df[fam_df["PID"] != "0"]
    n_has_father = has_father.shape[0]
    has_father_in_data = fam_df[(fam_df["PID"] != "0") & (fam_df["PID"].isin(fam_df["IID"]))]
    n_has_father_in_data = has_father_in_data.shape[0]
    n_missing_fathers=n_has_father-n_has_father_in_data
    n_kinship_clusters = fam_df.drop_duplicates(subset=["FID"]).shape[0]

    md_file.write(f"\n{n_kinship_clusters} kinship clusters")
    md_file.write(f"\n{n_has_mother} offspring with mother ID")
    md_file.write(f"\n{n_missing_mothers} mothers missing from dataset")
    md_file.write(f"\n{n_has_father} offspring with father ID")
    md_file.write(f"\n{n_missing_fathers} fathers missing from dataset")

    md_file.write(f"\n## Call rates")
    md_file.write(f"\n### Sample call rates")
    imiss_df = pd.read_csv(imiss, delim_whitespace=True)
    sample_call_rates = 1-imiss_df["F_MISS"]
    md_file.write(f"\nmin: {sample_call_rates.min()}")
    md_file.write(f"\nmax: {sample_call_rates.max()}")
    md_file.write(f"\nmedian: {sample_call_rates.median()}")
    outTrunk = mqc.plinkBase(output_filename)
    plt.figure()
    sample_call_rates.plot.hist(bins=30, alpha = 0.7, color='blue')
    plt.title("Sample call rates")
    plt.xlabel("Call rate")
    plt.ylabel("Counts")
    plt.savefig(outTrunk+'.sample_call_rates_histogram.png')
    md_image_syntax = f'![Call Rates Histogram]({outTrunk}.call_rates_histogram.png)\n'
    md_file.write(md_image_syntax)
    md_file.close()