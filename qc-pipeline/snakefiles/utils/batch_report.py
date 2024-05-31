import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mobaQcTools as mqc
import os
def write_report(output_filename, batch, module, file_trunk):
    md_file = open(output_filename, "a")
    md_file.write(f"# Batch report for batch {batch}, module {module}")
    fam_df = pd.read_csv(file_trunk + ".fam", delim_whitespace=True, header=None, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])
    md_file.write(f"\n## Samples overview")
    included_number_of_samples = fam_df.shape[0]
    md_file.write(f"\n{included_number_of_samples} samples")

    n_has_mother, n_has_mother_in_data, n_missing_mothers, n_mothers = parental_data("MID", fam_df)
    n_has_father, n_has_father_in_data, n_missing_fathers, n_fathers = parental_data("PID", fam_df)

    n_kinship_clusters = fam_df.drop_duplicates(subset=["FID"]).shape[0]

    md_file.write(f"\n<br>{n_kinship_clusters} kinship clusters")
    md_file.write(f"\n<br>{n_has_mother} offspring with mother ID")
    md_file.write(f"\n<br>{n_has_mother_in_data} offspring with mother in batch")
    md_file.write(f"\n<br>{n_mothers} mothers with offspring in batch")
    md_file.write(f"\n<br>{n_missing_mothers} mothers missing from batch")
    md_file.write(f"\n<br>{n_has_father} offspring with father ID")
    md_file.write(f"\n<br>{n_has_father_in_data} offspring with father in batch")
    md_file.write(f"\n<br>{n_fathers} fathers with offspring in batch")
    md_file.write(f"\n<br>{n_missing_fathers} fathers missing from batch")

    md_file.write(f"\n## Call rates")
    imiss_df = pd.read_csv(file_trunk + ".imiss", delim_whitespace=True)
    lmiss_df = pd.read_csv(file_trunk + ".lmiss", delim_whitespace=True)
    write_stats_and_histogram(md_file, "Sample call rates", 1-imiss_df["F_MISS"], output_filename, x_label="Call rate")
    write_stats_and_histogram(md_file, "SNP call rates", 1-lmiss_df["F_MISS"], output_filename, x_label="Call rate")

    md_file.write(f"\n## F_het")
    fhet_df = pd.read_csv(file_trunk + ".het", delim_whitespace=True)
    write_stats_and_histogram(md_file, "F_het", fhet_df["F"], output_filename, x_label="$F_{het}$", subheader= False)

    md_file.write(f"\n## Hardy-Weinberg P-values")
    hwe_df = pd.read_csv(file_trunk + ".hwe", delim_whitespace=True)
    write_stats_and_histogram(md_file, "Hardy-Weinberg P-values", hwe_df["P"], output_filename, x_label="P-value", subheader= False)

    md_file.write(f"\n## Sexcheck")
    sexcheck = pd.read_csv(file_trunk + ".sexcheck", delim_whitespace=True)
    ok_status = sexcheck[sexcheck["STATUS"] == "OK"]
    n_ok_status = ok_status.shape[0]
    md_file.write(f"\n{n_ok_status} out of {included_number_of_samples} OK<br>\n")
    write_sexcheck_table(md_file, sexcheck)
    
    md_file.write(f"\n### All samples")
    write_sexcheck_scatterplot(md_file, sexcheck, "all_F.png", os.path.dirname(output_filename), "F-statistics for all samples", groupby="SNPSEX")

    md_file.write(f"\n### PEDSEX Male")
    pedsex_male = sexcheck[sexcheck["PEDSEX"] == 1]
    write_sexcheck_scatterplot(md_file, pedsex_male, "male_F.png", os.path.dirname(output_filename), "F-statistics for PEDSEX Male")
    write_stats_and_histogram(md_file, "PEDSEX Male F-statistics", pedsex_male["F"], output_filename, x_label="F", subheader=True)
    
    md_file.write(f"\n### PEDSEX Female")
    pedsex_female = sexcheck[sexcheck["PEDSEX"] == 2]
    write_sexcheck_scatterplot(md_file, pedsex_female, "female_F.png", os.path.dirname(output_filename), "F-statistics for PEDSEX Female")
    write_stats_and_histogram(md_file, "PEDSEX Female F-statistics", pedsex_female["F"], output_filename, x_label="F", subheader=True)
    md_file.close()

def write_sexcheck_scatterplot(md_file, sexcheck, png_file, output_path, title, groupby = 'SNPSEX'):
    colors = {0: 'red', 1: 'green', 2: 'blue'}
    color_list = [colors[group] for group in sexcheck[groupby]]
    legend_handles = [mpatches.Patch(color=colors[1], label=f"{groupby} Male"), mpatches.Patch(color=colors[2], label=f"{groupby} Female"), mpatches.Patch(color=colors[0], label=f"{groupby} Unknown")]
    ax = sexcheck.plot.scatter("F", "YCOUNT", c=color_list, grid=True)
    ax.set_ylim([0, None])
    ax.legend(handles=legend_handles, loc='upper left')
    ax.set_title(title)
    plt.savefig(png_file)
    png_path = f"{output_path}/{png_file}"
    plt.savefig(png_path, dpi=200)
    md_image_syntax = f"<br><img src='{png_file}' width='700'/>"
    md_file.write(md_image_syntax)


def write_sexcheck_table(md_file, sexcheck):
    md_file.write("| PEDSEX | Total | SNPSEX Male | SNPSEX Female | SNPSEX Unknown | OK | Problem |\n")
    md_file.write("| ------ | ------ | ------ | ------ | ------ | ------ | ------ |\n")
    sexcheck_male = sexcheck[sexcheck["PEDSEX"] == 1]
    sexcheck_female = sexcheck[sexcheck["PEDSEX"] == 2]
    sexcheck_unknown = sexcheck[sexcheck["PEDSEX"] == 0]
    write_sexcheck_stats(md_file, sexcheck_male, "Male")
    write_sexcheck_stats(md_file, sexcheck_female, "Female")
    write_sexcheck_stats(md_file, sexcheck_unknown, "Unknown")


def write_sexcheck_stats(md_file, sexcheck, sex):
    n_total = sexcheck.shape[0]
    n_males = sexcheck[sexcheck["SNPSEX"] == 1].shape[0]
    n_females = sexcheck[sexcheck["SNPSEX"] == 2].shape[0]
    n_unknown = sexcheck[sexcheck["SNPSEX"] == 0].shape[0]
    n_ok = sexcheck[sexcheck["STATUS"] == "OK"].shape[0]
    n_problem = sexcheck[sexcheck["STATUS"] == "PROBLEM"].shape[0]
    md_file.write(f"| {sex} | {n_total} | {n_males} | {n_females} | {n_unknown} | {n_ok} | {n_problem} |\n")

def parental_data(parent_type, fam_df):
    has_parent=fam_df[fam_df[parent_type] != "0"]
    n_has_parent = has_parent.shape[0]
    has_parent_in_data = fam_df[(fam_df[parent_type] != "0") & (fam_df[parent_type].isin(fam_df["IID"]))]
    n_has_parent_in_data = has_parent_in_data.shape[0]
    n_missing_parents=n_has_parent-n_has_parent_in_data
    parents = fam_df[fam_df["IID"].isin(fam_df[parent_type])]
    n_parents = parents.shape[0]
    return n_has_parent, n_has_parent_in_data, n_missing_parents, n_parents

def write_stats_and_histogram(md_file, title, series, output_filename, x_label = "Value", subheader = True):
    if subheader:
        md_file.write(f"\n### {title}")
    md_file.write(f"\nmin: {series.min()}")
    md_file.write(f"\n<br>max: {series.max()}")
    md_file.write(f"\n<br>median: {series.median()}")
    plt.figure()
    series.plot.hist(bins=100, alpha = 0.7, color='blue')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel("Counts")
    path = os.path.dirname(output_filename)
    title_png = f'{title.replace(" ","_")}_histogram.png'
    title_path = f"{path}/{title_png}"
    plt.savefig(title_path, dpi=200)
    # md_image_syntax = f'\n<br>![]({title_png})'
    md_image_syntax = f"<br><img src='{title_png}' width='700'/>"
    md_file.write(md_image_syntax)