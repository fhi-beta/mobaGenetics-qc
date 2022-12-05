#!/usr/bin/python3

# 4.11.2019 This files is probably not in use. Kept for now as the logging done here is more general,
# but maybe not necessary

import subprocess
import logging

# Logging should only by initialised once ... Fix that later
logging.basicConfig(filename=mqc.globalLog,
                    filemode='a',
                    format='%(asctime)s :: %(name)s :: %(levelname)s :: %(lineno)d :: %(message)s',
                    level = logging.INFO)

logging.info("inputfiles " +  ' '.join(snakemake.input))
logging.info("to make output " +  ' '.join(snakemake.output))

# A little hack to reduce the number of params: From snakemake.input/output we know where the 
# Input and ouput-trunks used by plink
inTrunk =  mqc.plinkBase(snakemake.input[0])
outTrunk =  mqc.plinkBase(snakemake.output[0])

subprocess.run([snakemake.params.plink,
                "--bfile",inTrunk,
                "--update-ids", snakemake.input.recode_id ,
                "--out", outTrunk ,
                "--make-bed"
               ])


# Checking common patterns before and after plinkedits. Note that this could have been done with set
# intersections, but we are having fun with only absorbing one of the files in memory. You should call
# dictFromFile with the smallest file, as check_match will just interate on the file.
cols = [0,1]
editedSet = mqc.dictFromFile(outTrunk + ".fam", cols)       
matches = mqc.checkMatch(inTrunk + ".fam", editedSet, cols)

logging.warning('Results NOT logged to separate file ' + "Common lines: " + str(matches) + " based on " + str(len(editedSet)) + " Samples.  Columns checked  "+",".join(map(str,cols)))


