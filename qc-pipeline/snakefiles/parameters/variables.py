#!/usr/bin/python3

# These global variables can be shared only between Snakefiles
### snakemake_workflows initialization ########################################


# Executables
libdir = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'lib'))
bcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/bcftools-1.7/bcftools'))
vcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/vcftools-0.1.13/vcftools'))
hrc1000g = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir),'bin/HRC-1000G/HRC-1000G-check-bim.pl'))
beagle = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/beagle/beagle.01Mar24.d36.jar'))
conform_gt = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/beagle/conform-gt.24May16.cee.jar'))
plinklocal = Path(os.path.dirname(workflow.basedir))/'bin'/'plink-1.90b5.4'/'plink'
plink2local = Path(os.path.dirname(workflow.basedir))/'bin'/'plink2'/'plink2'
flashpca = Path(os.path.dirname(workflow.basedir))/'bin'/'flashpca_x86-64'
kinglocal = Path(os.path.dirname(workflow.basedir))/'bin'/'Linux-king'/'king'

# Resources
high_ld_regions_hg19 = Path(os.path.dirname(workflow.basedir))/'resources'/'high-ld-regions-hg19'
hrc_sites = Path(config['hrc_sites'])
mapfiles = Path(config['mapfiles'])
exclude_variants = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'resources/exclude_variants.txt')) # used to exclude variants with missing IDs

# Folder paths
tmp_path = Path(config['output_base']) / 'tmp'
runlog = Path(config['output_base']) / 'runlog.txt'
base = Path(config['output_base'])
hrc_ega = Path(config['hrc_ega'])
github_docs = Path(os.path.dirname(workflow.basedir)) / 'docs'
tmpMod1 = base/'mod1-data-preparation'
tmpMod2 = base/'mod2-genetic-relationship'
tmpMod3 = base/'mod3-population-clustering'
tmpMod4 = base/'mod4-good_markers'
tmpMod5 = base/'mod5-pre-phasing'
tmpMod6 = base/'mod6-imputation'
# tmpMod5 = base/'mod5-samples_unrelated'
# tmpMod6 = base/'mod6-phasing-preparation'
resultPath = base/'results'

### Batch settings
batches = ['snp001', 'snp002', 'snp003', 'snp007', 'snp008', 'snp009', 'snp010', 'snp011', 'snp012', 'snp014', 'snp015a', 'snp015b', 'snp016a', 'snp016b', 'snp017a', 'snp017b', 'snp017c', 'snp017d', 'snp017e', 'snp017f', 'snp018a', 'snp018b', 'snp018c', 'snp018d', 'snp018e']
batches = ['snp007', 'snp008', 'snp009', 'snp011', 'snp012', 'snp014', 'snp015a', 'snp015b', 'snp016a', 'snp016b', 'snp017a', 'snp017b', 'snp017c', 'snp017d', 'snp017e', 'snp017f', 'snp018a', 'snp018d'] # Checked up to module 3
batches = ['snp011', 'snp012', 'snp014', 'snp015a', 'snp015b', 'snp016a', 'snp016b', 'snp017a', 'snp017b', 'snp017c', 'snp017d', 'snp017e', 'snp017f', 'snp018a', 'snp018b']


### workflow settings ##################################
roles = ["founders", "offspring"]

