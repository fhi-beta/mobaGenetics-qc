#!/usr/bin/python3
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from plotnine import *

import re
import datetime
from pathlib import Path

def plotHist(dataFile,resultFile,column="name of the column"):
    """
    Very degenerated test
    """
    df = pd.read_csv(dataFile, sep="\t",
                 usecols = [column] )
    df = df.sort_values("10% GC")
    p = ggplot(data=df, mapping=aes(x='10% GC'))
    hist = p + geom_histogram()
    ggsave(plot=hist, filename=resultFile, dpi=600)
    return
        

def plinkBase(path):
    """
    plink often works on a trunk and creates extra files: a becomes a.fam, a.bed, a.bid
    This function reates this trunk, so /foo/bar/gazonk.fam becomes /foo/bar/gazonk
    The resulting trunk is typically used as input to plink
    Only the last '.' is removed so /foo.foo/bar/gazonk.x.fam will become /foo.foo/bar/gazonk.x
    """
    return  re.sub(r"\.\w*$","",path)


def lt(a,b):
    return a < b
def eq(a,b):
    return a==b
def gt(a,b):
    return a > b
def unknownComp(a,b):
    print ("Unknown compare operation")
    return 



# used to map symbols to actual comparision functions defined previously
compOperands = {
    "<": lt,
    "=": eq,
    ">": gt
}


def extractSampleList(innFile, sampleFile, subsetFile = "/dev/null", 
                      colName = "none", sep = '\t', condition = "<", treshold = 0, cols = {1}, ):
    """
    A typical preprocessor to utilities like plink. Efficient in the sense that it 
    is not reading the (huge) files into memory.
    Takes a csv file innFile (first line with headers, seperator is paramneter "sep") and produces 
    * a sample list on sampleFile (one columun, for now hardcoded to column 0)
    * a subset file with at set of numbered columns (cols) that also contains colName (see below)
    
    Only columns were the column colName matches the condition of treshhold will be written
    Returns (number of sample extracted, total samples)
    Restrictions: 
    * Assuming the first column in innFile is samplenumber - Only the first column will be written to sampleFile
    * Assumes only one match for the regExp colName
    * Columns are ordered in the same way as the innFile : {0,1} is identical to {1,0,0,1}
    * Outputfiles use the same separator as used for innFile
    """ 
    with open(innFile) as fp:
        
        # Identifying header column 
        line = fp.readline()
        allcols = line.split(sep)
        regex = re.compile(colName)
        indx = [i for i, item in enumerate(allcols) if re.search(regex, item)]
        # print(f"{colName} is column {indx}")
        matches = 0
        compare = compOperands.get(condition, unknownComp)
        if compare == unknownComp :
            print (f"Cannot compare using: '{condition}'")
            return
        # These lines are data, indx[0] is the index  of the column we found interesting
        sample = open(sampleFile,"w+")
        subset = open(subsetFile,"w+")
        for line in enumerate(fp):
            allcols = line[1].split(sep)
            if compare(float(allcols[indx[0]]) , treshold): 
                # File with only sample number
                sample.write(allcols[0] + "\n")
                # File with more columns, and always include the treshold column
                # subsetc is a relevant subset of the columns 
                subsetc = [allcols[index] for index in cols]                
                subset.write(sep.join(map(str,subsetc)) + sep + allcols[indx[0]] + "\n" )
                matches += 1
    totalLines = line[0] + 1 # enumerates from 0, and we read a line manually
    sample.close()
    subset.close()
    return (matches,totalLines)

    

    
def dictFromFile(fil,cols=[0,1]):
    """
    Creates a dictionary from the concatenation of the values in the columns found in fil.
    fil is expected to be a csv file with whitespace as delimiters
    First Column is number 0.
    (Hint for when this is used to compare the columns of two files: You can reduce memory usage by 
    making a dictionary of the smallest file and searchin/iterating for matches through the largest

    Will print error/warning if dictionary is empty or if replicates found
    """
    # pandas is overkill here, but since we will use it anyway ...
    # Grab the relevant columns only
    df = pd.read_csv(fil, usecols=cols, delim_whitespace=True, header=None)
    # concat them as strings, handling NA
    df = pd.Series(df.fillna('').values.tolist()).map(lambda x: ''.join(map(str,x)))
    # ... and finally make a dictionary. Note that we will here also count identical lines.
    countDict = dict()
    max = 0
    for i in df:  
        countDict[i] = countDict.get(i,0) + 1
        if ( countDict[i]) > max : 
            max = countDict[i]
            sample = i
    if len(countDict) == 0 : print("ERROR: 0 sized dictionary after reading ",fil)
    if max>1 : print("WARNING: sample " + str(sample) + " found " + str(max) + "times. Might not be the only nonunique")
    return countDict

def checkMatch(fil,dic,cols=[0,1]):
    """
      For each line in file  fil, extract columns cols and check if their concatenation exist in dictionary dic
      Returns (number of matches , number of lines checked)
      File must be a csv file with whitespace as delimiters
    """
    matches = 0
    lines = 0
    for line in open(fil): 
        lines += 1
        allcols = line.split()
        subsetc = [allcols[index] for index in cols]
        # concatenating so we can look up the strings in the existing dictionary
        key = "".join(map(str,subsetc))
        if (dic.get(key,0) > 0): matches += 1
    return (matches,lines)


def countCsvDiff(bigfile, smallfile, cols = [0,1]):
    """
    Every now and then we have two csv files that need to be compared, but only certain columns (default the two first)
    This does exactly that, and counts the number of matches
    Most efficient if you pass the largest file's name in bigfile
    Returns the (matches, number_of_lines_in_big_file - 1)
    Note that if bigfile has a header, this is the number of datarows in the a csv file
    """
    smallDict = dictFromFile(smallfile, cols)       
    matches = checkMatch(bigfile, smallDict, cols)

    # logging.warning('Results NOT logged to separate file ' + "Common lines: " + str(matches) + " based on " + str(len(smallDict)) + " Samples.  Columns checked  "+",".join(map(str,cols)))
    return matches

def log(logfile, message = "Nothing logged", mode = "a" ):
    """
    Placeholder until we figure out how logging really needs to be done. Default appends message to logfile
    """
    with open(logfile, mode) as myfile:
        myfile.write('{:%Y-%m-%d %H:%M:%S} '.format(datetime.datetime.now()) + message)
    return

def resultLog(logfile, message = "Nothing logged", mode = "w+" ):
    """
    Placeholder until we figure out how results will be reporte. Default creates a new message to logfile
    No formating is done, in contrary to log()
    Currently it logs to logfile, but also to a global file (always in append-mode) where it also shows the path of the result.
    The (crude) idea of the global file is to make it easier to make a comple (and prettified) report. The name should be a variable somehow
    Have not tested this with threads, it will probably fail unless you reset the file and use append mode
    """
    with open(logfile, mode) as myfile:
        myfile.write(message)

    # Dirty - hardcoded name here ...
    mergedResults = Path(logfile).parent/"mergedResults.txt"
    with open(mergedResults,"a") as myfile:
        myfile.write( f"(ref {logfile} ) " + message)
    return

def _make_gen(reader):
    """
    Used to linecount buffered and fast
    """ 
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)


def lineCount(filename):
    """
    Efficient linecounting by buffered reads
    """
    f = open(filename, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum( buf.count(b'\n') for buf in f_gen )
