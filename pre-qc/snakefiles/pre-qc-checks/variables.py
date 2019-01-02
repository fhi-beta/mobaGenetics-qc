### snakemake_workflows initialization ########################################
libdir = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'lib'))
#plinklocal = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'software/plink-1.90b5.3/plink'))


### workflow settings ##################################
chrom = list(range(1,23))


### generate paths ###################################
if not os.path.exists(config['output_base']):
    os.makedirs(config['output_base'])

if not os.path.exists(config['tmp_path']):
    os.makedirs(config['tmp_path'])

