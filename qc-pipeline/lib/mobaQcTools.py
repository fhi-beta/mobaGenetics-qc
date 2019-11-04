#!/usr/bin/python3
import pandas as pd
import re
import datetime

def plinkBase(path):
    """
    plink often works on a trunk and creates extra files: a becomes a.fam, a.bed, a.bid
    This function reates this trunk, so /foo/bar/gazonk.fam becomes /foo/bar/gazonk
    The resulting trunk is typically used as input to plink
    Only the last '.' is removed so /foo.foo/bar/gazonk.x.fam will become /foo.foo/bar/gazonk.x
    """
    return  re.sub(r"\.\w*$","",path)

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
    if max>1 : print("WARNING: sample " + str(sample) + " found " + str(max) + "times. Mignt not be the only nonunique")
    return countDict

def checkMatch(fil,dic,cols=[0,1]):
    """
      For each line in file  fil, extract columns cols and check if their concatenation exist in dictionary dic
      Returns the number of matches
      File must be a csv file with whitespace as delimiters
    """
    matches = 0
    for line in open(fil): 
        allcols = line.split()
        subsetc = [allcols[index] for index in cols]
        # concatenating so we can look up the strings in the existing dictionary
        key = "".join(map(str,subsetc))
        if (dic.get(key,0) > 0): matches += 1
    return matches


def countCsvDiff(bigfile, smallfile, cols = [0,1]):
    """
    Every now and then we have two csv files that need to be compared, but only certain columns (default the two first)
    This does exactly that, and counts the number of matches
    Most efficient if you pass the largest file's name in bigfile
    """
    smallDict = dictFromFile(smallfile, cols)       
    matches = checkMatch(bigfile, smallDict, cols)

    # logging.warning('Results NOT logged to separate file ' + "Common lines: " + str(matches) + " based on " + str(len(smallDict)) + " Samples.  Columns checked  "+",".join(map(str,cols)))
    return matches

def log(logfile, message):
    """
    Placeholder until we figure out how logging really needs to be done. Appends message to logfile
    """
    with open(logfile,"a") as myfile:
        myfile.write('{:%Y-%m-%d %H:%M:%S} '.format(datetime.datetime.now()) + message)
    return

