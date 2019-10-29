#!/usr/bin/python3

### snakemake_workflows initialization ########################################
libdir = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'lib'))
bcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/bcftools-1.7/bcftools'))
vcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/vcftools-0.1.13/vcftools'))
plinklocal = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/plink-1.90b5.4/plink'))
tmp_path = os.path.join(config['output_base'],'tmp')
runlog = config['output_base'] + '/runlog.txt'
base = config['output_base']

### workflow settings ##################################
# chrom gen. list 1-23 (24 not in list)
chrom = list(range(1,24))
chromx = list(range(1,23)) + ['X']
ROLES = ["founders", "offspring"]
#mod1 = os.path.join(config['output_base'], 'mod1-data-preparation')
#modx = os.path.join(config['output_base'], 'modx-shaping-preparation')

### generate paths ###################################
if not os.path.exists(config['output_base']):
    os.makedirs(config['output_base'])

if not os.path.exists(tmp_path):
    os.makedirs(tmp_path)

#for p in [mod1,modx]:
#	if not os.path.exists(p):
#		os.makedirs(p)



