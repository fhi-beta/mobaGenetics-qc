#!/usr/bin/python3

# These global variables can be shared only between Snakefiles
### snakemake_workflows initialization ########################################


# Executables
libdir = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'lib'))
bcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/bcftools-1.7/bcftools'))
vcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/vcftools-0.1.13/vcftools'))
hrc1000g = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir),'bin/HRC-1000G/HRC-1000G-check-bim.pl'))
beagle = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/beagle/beagle.06Aug24.a91.jar'))
conform_gt = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/beagle/conform-gt.24May16.cee.jar'))
plinklocal = Path(os.path.dirname(workflow.basedir))/'bin'/'plink-1.90b5.4'/'plink'
plink2local = Path(os.path.dirname(workflow.basedir))/'bin'/'plink2'/'plink2'
flashpca = Path(os.path.dirname(workflow.basedir))/'bin'/'flashpca_x86-64'
kinglocal = Path(os.path.dirname(workflow.basedir))/'bin'/'Linux-king'/'king'

# Resources
resources_folder = Path(os.path.dirname(workflow.basedir))/'resources'
high_ld_regions_hg19 = resources_folder/'high-ld-regions-hg19.txt'
haploid_x_file = resources_folder/'haploid_x.txt'
hrc_sites = Path(config['hrc_sites'])
sites_1000g = Path(config['1000g_sites'])
mapfiles = Path(config['mapfiles'])
exclude_variants = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'resources/exclude_variants.txt')) # used to exclude variants with missing IDs
blacklisted_variants = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'resources/variant_blacklist.txt'))

# Folder paths
tmp_path = Path(config['output_base']) / config['release'] / 'tmp'
runlog = Path(config['output_base']) / config['release'] / 'runlog.txt'
base = Path(config['output_base']) / config['release']
archive_base = Path(config['archive_output_base']) / config['release']
archive3_base = Path(config['archive3_output_base']) / config['release']
archive3_base_2 = Path(config['archive3_output_base']) / config['release2']
hrc_ega = Path(config['hrc_ega'])
hrc_bref = Path(config['hrc_bref'])
hrc_vcf = Path(config['hrc_vcf'])
github_docs = Path(os.path.dirname(workflow.basedir)) / 'docs' / config['release']

github_docs_2 = Path(os.path.dirname(workflow.basedir)) / 'docs' / config['release2']
tmpMod1 = base/'mod1-data-preparation'
tmpMod2 = base/'mod2-genetic-relationship'
tmpMod3 = base/'mod3-population-clustering'
tmpMod4 = base/'mod4-good_markers'
tmpMod5 = base/'mod5-pre-phasing'
tmpMod6 = base/'mod6-imputation'
tmpMod7 = base/'mod7-post-imputation'
tmpMod8 = archive3_base/'mod8-release_annotation'
tmpMod8_2 = archive3_base_2/'mod8-release_annotation'
#tmpMergeTest = base/'merge_test'
# use archive for these for now, since /work is full:
#tmpMod6Archive = archive_base/'mod6-imputation'
#tmpMod7 = base/'mod7-post-imputation'
#tmpMod7 = Path(config['output_base']) / '2024.09.23' / 'mod7-post-imputation'
# tmpMod72 = Path(config['output_base']) / '2024.09.23' / 'mod7-post-imputation'
#tmpMod8 = Path(config['output_base']) / config['release2'] / 'mod8-release_annotation'
release_folder = Path(config['release_base']) / config['release'] # Path(config['release_base']) / config['release']
release_folder_2 = Path(config['release_base']) / config['release2'] # Path(config['release_base']) / config['release']
release_base_name = "moba_genotypes_" + config['release'] # "moba_genotypes_" + config['release']

n_samples = config['n_samples']
# tmpMod5 = base/'mod5-samples_unrelated'
# tmpMod6 = base/'mod6-phasing-preparation'
resultPath = base/'results'

### Batch settings
# batches = ['snp001', 'snp002', 'snp003', 'snp007', 'snp008', 'snp009', 'snp010', 'snp011', 'snp012', 'snp014', 'snp015a', 'snp015b', 'snp016a', 'snp016b', 'snp017a', 'snp017b', 'snp017c', 'snp017d', 'snp017e', 'snp017f', 'snp018a', 'snp018b', 'snp018c', 'snp018de']
# Notes:
# - Batches 18ab seem to be parent/child only, check the selection criteria.
# - Batch 13 lacks documentation and needs to be checked for overlap with other batches before being added. Gutorm notes that it probably cannot be included.
# - Batch 19 needs to be checked for overlap with other batches before being added. Gutorm notes that it probably does not need to be included.
# batches_debug = ['snp001', 'snp002', 'snp003', 'snp007', 'snp008', 'snp009', 'snp010', 'snp011', 'snp012', 'snp014', 'snp015a', 'snp015b', 'snp016a', 'snp016b', 'snp017a', 'snp017b', 'snp017c', 'snp017d', 'snp017e', 'snp017f', 'snp018a', 'snp018b', 'snp018c', 'snp018de']
batches =  config['batches']
chrs = config['chrs']
# chrs_debug = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16']
# chrs_debug = [str(i) for i in range(1, 23)]