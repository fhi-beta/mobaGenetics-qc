#!/usr/bin/python3
import pandas as pd

globalLog = "./morn"

def dictFromFile(fil,cols=[0,1]):
    """
    Creates a dictionary from the concatenation of the values in the columns found in fil.
    fil is expected to be a csv file with whitespace as delimiters
    First Column is number 0.
    (Hint for when this is used to compare the columns of two files: You can reduce memory usage by 
    making a dictionary of the smallest file and searchin/iterating for matches through the largest

    Will print error/warning if dictionary is empty or if replicates found
    """
    # print("**** dictFromFile "+ fil)
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
    #print("**** checkMatch "+ fil)
    matches = 0
    for line in open(fil): 
        allcols = line.split()
        subsetc = [allcols[index] for index in cols]
        # concatenating so we can look up the strings in the existing dictionary
        key = "".join(map(str,subsetc))
        if (dic.get(key,0) > 0): matches += 1
    return matches



# cols = [0,1,2]

# # Checking common patterns before and after plinkedits. Note that this could have been done with set
# # intersections, but we are having fun with only absorbing one of the files in memory. You should call
# # dictFromFile with the smallest file, as check_match will just interate on the file.

# orig = dictFromFile("foo.fam",cols)
# matches = check_match("bar.fam",orig,cols)

# print("Common lines: " + str(matches) + " Columns checked")
# print(cols)



# # print (dictFromFile.__doc__)

