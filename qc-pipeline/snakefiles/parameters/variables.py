#!/usr/bin/python3
from pathlib import Path
from shutil import copyfile
# numerical & plots
import pandas as pd
import numpy as np
import plotnine as p9
import matplotlib
import math
import sys
import yaml
import io
import shutil
import os
matplotlib.use('Agg')  

# These global variables can be shared only between Snakefiles 
### snakemake_workflows initialization ########################################

libdir = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'lib'))
bcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/bcftools-1.7/bcftools'))
vcftools = os.path.abspath(os.path.join(os.path.dirname(workflow.basedir), 'bin/vcftools-0.1.13/vcftools'))

# Local binaries
plinklocal = Path(os.path.dirname(workflow.basedir))/'bin'/'plink-1.90b5.4'/'plink'
plink2local = Path(os.path.dirname(workflow.basedir))/'bin'/'plink-1.90b5.4'/'plink'
flashpca = Path(os.path.dirname(workflow.basedir))/'bin'/'flashpca_x86-64'
kinglocal = Path(os.path.dirname(workflow.basedir))/'bin'/'Linux-king'/'king'

# Folder paths
tmp_path = Path(config['output_base']) / 'tmp'
runlog = Path(config['output_base']) / 'runlog.txt'
base = Path(config['output_base'])
github_docs = Path(os.path.dirname(workflow.basedir)) / 'docs'
tmpMod1 = base/'mod1-data-preparation/'
tmpMod2 = base/'mod2-data-preparation/'
tmpMod3 = base/'mod3-good-markers/'
# tmpMod2 = base/'mod2-genetic-relationship/'
# tmpMod3 = base/'mod3-population-clustering/'
tmpMod4 = base/'mod4-samples-unrelated/'
tmpMod5 = base/'mod5-shaping-preparation/'
resultPath = base/'results/'

### Batch settings
# batches = ['snp001', 'snp002', 'snp003', 'snp007', 'snp008', 'snp009', 'snp010', 'snp011', 'snp012', 'snp014', 'snp015a', 'snp015b', 'snp016a', 'snp016b', 'snp017a', 'snp017b', 'snp017c', 'snp017d', 'snp017e', 'snp017f', 'snp018a', 'snp018b', 'snp018c', 'snp018d', 'snp018e']
batches = ['snp012', 'snp014']
batches = ['snp014']


### workflow settings ##################################
chrom = list(range(1,24))
chromx = list(range(1,23)) + ['X']
ROLES = ["founders", "offspring"]

