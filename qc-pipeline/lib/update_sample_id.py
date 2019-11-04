#!/usr/bin/python3
import subprocess
import logging
import mobaQcTools as mqc  # local moba package
# Logging should only by initialised once ... Fix that later
logging.basicConfig(filename=mqc.globalLog,
                    filemode='a',
                    format='%(asctime)s :: %(name)s :: %(levelname)s :: %(lineno)d :: %(message)s',
                    level = logging.INFO)

logging.info("inputfiles " +  ' '.join(snakemake.input))
logging.info("to make output " +  ' '.join(snakemake.output))

# For now - not checking result yet. Might be wrapptd later
subprocess.run([snakemake.params.plink,
                "--bfile",snakemake.params.in_bedset,
                "--update-ids", snakemake.input.recode_id ,
                "--out", snakemake.params.out_bedset ,
                "--make-bed"
               ])


# Checking common patterns before and after plinkedits. Note that this could have been done with set
# intersections, but we are having fun with only absorbing one of the files in memory. You should call
# dictFromFile with the smallest file, as check_match will just interate on the file.
cols = [0,1]
editedSet = mqc.dictFromFile(snakemake.params.in_bedset+".fam", cols)
matches = mqc.checkMatch(snakemake.params.out_bedset+".fam", editedSet, cols)

logging.warning("WTF")
logging.warning('Results NOT logged to separate file ' + "Common lines: " + str(matches) + " based on " + str(len(editedSet)) + " Samples.  Columns checked  "+",".join(map(str,cols)))
subprocess.run(["touch","foo"])


