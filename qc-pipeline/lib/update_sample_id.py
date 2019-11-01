#!/usr/bin/python3
import os
import logging
import mobaQcTools as mqc  # local moba package

print ("at least i tried")

# Logging should only by initialised once ... Fix that later
logging.basicConfig(filename=mqc.globalLog,
                    filemode='a',
                    format='%(asctime)s :: %(name)s :: %(levelname)s :: %(lineno)d :: %(message)s',
                    level = logging.INFO)

print ("at least i tried again")
# logging.warning('Nothing to warn about')
logging.info("inputfiles " +  ' '.join(snakemake.input))
logging.info("to make output " +  ' '.join(snakemake.output))

plinkcmd = snakemake.params.plink \
               + " --bfile " + snakemake.params.in_bedset \
               + " --update-ids " + snakemake.input.recode_id \
               + " --out " + snakemake.params.out_bedset \
               + " --make-bed"
# use subpreocess
#os.system(plinkcmd)
print ("Would run " + plinkcmd)


# Checking common patterns before and after plinkedits. Note that this could have been done with set
# intersections, but we are having fun with only absorbing one of the files in memory. You should call
# dictFromFile with the smallest file, as check_match will just interate on the file.
cols = [0,1]
editedSet = mqc.dictFromFile(snakemake.params.in_bedset+".fam", cols)
matches = mqc.checkMatch(snakemake.params.out_bedset+".fam", editedSet, cols)

print("Common lines: " + str(matches) + " based on " + str(len(editedSet)) + " Samples.  Columns checked  ")
print(cols)




