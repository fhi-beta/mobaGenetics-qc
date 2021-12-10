Based on prior QC work by Christian M. Page and Haakon E. Nustad and
the bhklab/IlluminaEPICmethylation pipeline on github by Christopher Eeles and Joanna Pryzbyl.
Developed by Hanna Zdanowicz, Haakon E. Nustad and Gutorm Thomas Høgåsen.

General set-up steps (outside of TSD):
Clone the git repository. every command must be run from inside the main directory.
use python (version 3) to create virtual environment, and install snakemake.

if you do not want to use renv, make sure renv/, renv.lock and .Rprofile are located in Setup folder and change config.yaml renv variable to FALSE


Set up on TSD:

Unpack the zipped pipeline folder at prefered destination on TSD.
move within folder and use python3 to set up virtual environment and install snakemake:
$python -m venv .venv
$pip install snakemake

move renv/ folder, renv.lock and .Rprofile out of Setup folder to root folder

change config.yaml renv variable to TRUE

for installation of packages into local R library:
Start interative session with R. This will trigger renv, which is version control and package management software for R.

# Start interactive R:
R

#trigger dependency installation if init() did not:
renv::restore()

# exit R
q()

Create some folders relevant for the pipeline. 
# Create folders in root:
mkdir {plots,logs,qc_results,processed_data} 

Go through config.yaml and change parameters to your specific file paths and names. Especially data directory (idat and sample sheet file path), data type (450k or EPIC). Additionally, add file in doc folder: manual_remove_samples_bad_density.txt. good to go!


Pre-quality control:
important pre-step: the sample sheet needs the beadchip array column to have column name: Sentrix_ID, 
and the position on the array to have column name: Sentrix_Position
This is important for minfi to automatically recognize correct columns to identify the .idat files.

All samples you wish to run QC for must be listed in your sample_sheet.csv, which must be located in the folder you specify in config.yaml. It should be placed together with the .idat files. This .csv file is recognized automatically with the minfi function applied to read in the sample_sheet. A current step implemented detects if you have a mix of children and older individuals, which will create an error with an error message, encouraging the user to split his/hers samples into two .csv files, and running these seperately. This is because the quality of the QC will increase if such biologically related differences are not present.


How to use:

$ source .venv/bin/activate
This activates the virtual environment, enabling snakemake

(if you wish to jump out of the virtual enviroment, do:
$ deactivate
)

the pipeline is run with the following commands:
$ snakemake --core 1 first_part
$ snakemake --core 1 second_part

it is organized in 3 parts, where one should inspect the plots produced located in plots/ after part 1 and 2. 
Based on the plots produced in first and second part, you should add sample ids to the file input/manual_remove_samples_bad_density.txt if you think they should be removed. Remove samples showing especially outlying patterns compared to the rest, based on the boxplot of control_probes and genome wide density plot.

After this, the third part can be run:
$ snakemake --core 1 third_part

At the end, a clean up rule can be run to remove unnecessary data:
$ snakemake --core 1 clean_up


The end matrices with data:
processed_data/:
{analysis_name}_12_Noob_beta_matrix_sex_chr.csv	(beta values from CpGs at sex chromosomes; rows - CpGs, columns - samples)
{analysis_name}_12_Noob_beta_matrix_sex_chr.rds	(same as above, but stored as .rds file)
{analysis_name}_13_beta_values.csv			(beta values from autosomal CpGs before BMIQ; rows - CpGs, columns - samples)
{analysis_name}_13_beta_values.rds			(same as above, but stored as .rds file)
{analysis_name}_13_bmiqed_beta_values.csv		(beta values from autosomal CpGs after BMIQ; rows - CpGs, columns - samples)
{analysis_name}_13_bmiqed_beta_values.rds		(same as above, but stored as .rds file)

these are stored as .csv files. To load into R, do the following command in R:
data = read.table("{analysis_name}_12_Noob_beta_matrix_sex_chr.csv", sep = ",", header = TRUE, row.names = 1, check.names = FALSE)

additional files with valueable information:
logs/:
The Rout files are log files containing information about the R session and the functions applied to the data. They contain information of what is done to the data during the different scripts. Additionally, if an error occurs in a script, the error message is given in the associated log file. These files are primarily meant for debugging of the pipeline, but can also be inspected to understand the order of when the different functions and QC steps have been applied.

plots/:
{analysis_name}_3_boxplot_control_probes.pdf (boxplot of log2()-values from intensity measures from control probes for each sample)
{analysis_name}_3_control_probe_PCA_plot.pdf (PC2 against PC1 from PCA of log2()-values from intensity measures from control probe)
{analysis_name}_6_rgSet_vs_Noob_density_plots.pdf (genomewide density from RAW methylation values compared to after ssNoob is applied. Estimated from a random (set.seed) selection of 50000 CpGs.)
{analysis_name}_8_qc_plot.pdf (minfi qc plot - median intensity from methylated vs unmethylated intensities. A mean of these being less than 10.5 indicates bad sample/low signal)
{analysis_name}_14_BMIQ_density_plots.pdf (similar to {analysis_name}_6_rgSet_vs_Noob_density_plots.pdf, but with genomewide density after BMIQ is applied as well)
{analysis_name}_15_histogram_of_differences.pdf (Comparison of RAW vs ssNoob methylation values. The difference of the two matrices are calculated, then the rowMeans (plot 1) and rowSds (plot 2) are calculated and plotted as histograms. The same is done between ssNoob methylation values and BMIQ methylation values. These plots are for documenting the difference after each major normalization method is applied.)

qc_results/:
{analysis_name}_2_probes_removed.csv (2 column matrix of CpG-ids removed and why: "cross_hybridizing", "polymorphic" (SNP influenced), "high_detection_p")
{analysis_name}_2_SNP_Betas.rds (matrix of beta values for the 65 (450k) or 60 (EPIC) SNPs on the array. rows - rs id, column - samples)
{analysis_name}_3_control_probe_PCA.rds (The PCA object of PCA analysis from control probes)
{analysis_name}_4_bisulphite_conversions.csv (The calculated bisulphite conversion rate for each sample)
{analysis_name}_4_num_probes_with_proportion_failed_samples_p0.01.csv (for each sample, the fraction of failed probes and the total of failed probes, compared to the detection p value cut off (0.01))
{analysis_name}_4_probes_failed_per_sample_p0.01.csv (for 5 different proportions (0.5, 0.1, 0.05, 0.01, 0.0033), how many CpGs have a bad detection p value for more than the give proportion of samples)
{analysis_name}_8_filtered_samples.csv (2 column matrix of sample ids removed and why: "NA_control_probes" (too many NA/INF values), "bisulphite" (low bisulphite conversion), "low_median_meth_or_unmeth_channel" (less than 10.5 mean median meth and unmeth log2-intensity), "manual_removed_bad_density" (removed due to having outlying values in either the control_probe_boxplot or genome-wide density), "pvalue" (removed samples with less than 90% detected probes based on detection pvalue cutoff))
{analysis_name}_9_estimated_cell_proportions.csv (matrix with samples as rows and cell type proportions estimated as columns)
{analysis_name}_11_added_sex_prediction_to_pheno.csv (sample sheet data with additional columns for predicted sex, median intensity across Y chromosome probes, and median intenisty across X chromosome probes)


Quality Control Pipeline:
- load iDat files to RGset (all information in minfi object format)

- filter probes from RGset
	- poor quality probes (detection p-value larger than 0.01 for more than 5 % of the samples)
	- cross-hybridizing probes (probes shown to have more than 1 possible binding site in the genome)
	- non-CpG methylation (CpH, (H = A, C or T))
	- probes influenced by SNPs with higher than 0.05 MAF. SNP placed in interrogated CpG or at the neighboring positions

- create boxplot and PCA plot of control probes:
	- PCA pot: PC2 against PC1 from PCA of approximately 1200 Redn and Green channel control probes
		- meant to use as documentation and to get intuition of possible amount of outliers
	- Boxplot of control probe values for each sample
		- if a sample has outlying pattern, such as very narrow range of values -> remove sample

- filter samples with:
	- more than 25% missing/infinite values for control probes
	- more than 10% probes with high detection p value (> 0.01)
	- less than 90% bisulphite conversion rate (calculated with bscon() function, which uses bisulphite conversion control probes to estimate conversion rate)
	
- preprocess with ss-Noob:
	- single-sample-(normal-exponential out-of-band) method used for background subtraction and dye-bias correction. It uses out-of-band probes to estimate a normal-exponential background signal for each color channel, and normalizes the signal from regular probes based on these estimated distributions to obtain a purer signal where the background is removed. This algorithm only uses information within each sample to correct each sample. Has been shown to be one of the better background subtraction and dye-bias correction methods: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5408810/
	- from the authors of ssNoob: "We note that on the Beta value scale, there is no difference between values returned by Noob or ssNoob. Differences are confined to the methylated and unmethylated signals."

- plotting of genomewide densities og RAW methylation data and after ssNoob is applied
	- if a sample has clear outlying distribution from the others, remove sample

- removing samples based on outlying pattern

- check median methylated intensity against median unmethylated intensity, and remove samples with (median_M + median_U)/2 < 10.5. Indicates low signal from either one or both channels.

- cell type proportion estimation. Uses either Blood reference set or CordBlood reference set, dependening on input to config.yaml.
Blood reference set: DOI: 10.1371/journal.pone.0041361
CordBlood reference set: doi: 10.1080/15592294.2016.1161875

- sex prediction and removement of sex discordant samples between given sex and predicted sex
	- measures the median intensity values from combined methylated and unmethylated signal from Y chromsome probes and the same for X chromosome probes. If the differences between these medians is more than 2 (less than -2) a female is predicted.
	
- seperate out sex chromosome measurements from autosomal measurements

- do BMIQ: probe type normalisation. the Beadchips have 2 different probe designs, I and II. These have a clear difference in dynamic range for the obtained methylation measurements. BMIQ projects the global distribution of probe type II measurements onto the distribution of probe type I, meaning that the methylation values from probe type II measuements are slightly shifted. This is especially a usual step when analysing regions of data, i.e. looking for differentially methylation regions. It is also believed that the dynamic range of probe type I is more correct than probe type II. BMIQ also calculates only within each sample.

- plotting og genomewide distribution of methylation values after BMIQ is applied.

- plot histogram of differences between raw methylation values and after ssNoob is applied, and between methylation values after ssNoob is applied and BMIQ is applied. These plots are for documenting the difference after each major normalization method is applied.

- clean up and delete unnecassry data stored at intermediate steps


More information:

detection p value:
a detection p value greater than 0.01 should not be trusted. calculated by comparing the total dna signal (methylated + unmethylated) for each position and sample to the background signal level. The background is estimated using the negative control probes, assuming a normal distribution. calculations are performed on the original (non-log) scale.

negative control probes: specifically design not to matchthe human genome, reads/signal from these probes in red and green channel to estimate the background signal distribution.

out of band probes: for a given type I bead, the intensity from the unused color channel has been proposed as a means of estimating background signal, and termed the out-of-band intensity. 



