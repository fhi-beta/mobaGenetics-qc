#!/usr/bin/python3
# import shutil
import os
import logging


print("morn")
logging.basicConfig(filename=snakemake.log[0],
                    filemode='a',
                    format='%(asctime)s :: %(name)s :: %(levelname)s :: %(lineno)d :: %(message)s',
                    level = logging.DEBUG)

logging.warning('Nothing to warn about')
logging.info("inputfiles " +  ' '.join(snakemake.input))
logging.info("to make output " +  ' '.join(snakemake.output))

plinklocal = snakemake.params.plink

#plinklocal="ls"

print ('Doing nothing ' + plinklocal)
os.system(plinklocal)

#shutil.copy(snakemake.input[0], snakemake.output[0])

#         $plinklocal \
#             --bfile {params.in_bedset} \
#             --update-ids {input.recode_id} \
#             --make-bed \
#             --out {params.tmp_path}/rawbed-id-updated


