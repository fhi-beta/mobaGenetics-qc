import pandas as pd
import matplotlib.pyplot as plt
import mobaQcTools as mqc
import os
def write_report(output_filename, batch, file_trunk):
    md_file = open(output_filename, "a")
    md_file.write(f"# Pre-imputation report for batch {batch}")
    fam_df = pd.read_csv(file_trunk + ".fam", delim_whitespace=True, header=None, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])
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

    md_file.write(f"\n<br>{n_kinship_clusters} kinship clusters")
    md_file.write(f"\n<br>{n_has_mother} offspring with mother ID")
    md_file.write(f"\n<br>{n_missing_mothers} mothers missing from dataset")
    md_file.write(f"\n<br>{n_has_father} offspring with father ID")
    md_file.write(f"\n<br>{n_missing_fathers} fathers missing from dataset")

    md_file.write(f"\n## Call rates")
    imiss_df = pd.read_csv(file_trunk + ".imiss", delim_whitespace=True)
    lmiss_df = pd.read_csv(file_trunk + ".lmiss", delim_whitespace=True)
    write_stats_and_plot(md_file, "Sample call rates", 1-imiss_df["F_MISS"], output_filename, x_label="Call rate")
    write_stats_and_plot(md_file, "SNP call rates", 1-lmiss_df["F_MISS"], output_filename, x_label="Call rate")

    md_file.write(f"\n## F_het")
    fhet_df = pd.read_csv(file_trunk + ".het", delim_whitespace=True)
    write_stats_and_plot(md_file, "F_het", fhet_df["F"], output_filename, x_label="$F_{het}$", subheader= False)

    md_file.write(f"\n## Hardy-Weinberg P-values")
    hwe_df = pd.read_csv(file_trunk + ".hwe", delim_whitespace=True)
    write_stats_and_plot(md_file, "Hardy-Weinberg P-values", hwe_df["P"], output_filename, x_label="P-value", subheader= False)

    md_file.close()

def write_stats_and_plot(md_file, title, series, output_filename, x_label = "Value", subheader = True):
    if subheader:
        md_file.write(f"\n### {title}")
    md_file.write(f"\nmin: {series.min()}")
    md_file.write(f"\n<br>max: {series.max()}")
    md_file.write(f"\n<br>median: {series.median()}")
    plt.figure()
    series.plot.hist(bins=30, alpha = 0.7, color='blue')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel("Counts")
    path = os.path.dirname(output_filename)
    title_png = f'{title.replace(" ","_")}_histogram.png'
    title_path = f"{path}/{title_png}"
    plt.savefig(title_path, dpi=200)
    md_image_syntax = f'\n<br>![]({title_png})'
    md_file.write(md_image_syntax)