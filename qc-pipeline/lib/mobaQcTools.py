#!/usr/bin/python3
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from plotnine import *
# While plotting in a non-interactive environment (xterm ...) we need this:
import matplotlib
matplotlib.use('Agg')  

import yaml
import json
import re
import datetime
from pathlib import Path
import subprocess
import inspect    # to find stack-info, such as the functions name
import os
from datetime import datetime
from shutil import copyfile

# This is kinda ugly
cwd = Path.cwd()
snake_dir = cwd/"../snakefiles"
# snake_dir = Path("/mnt/work/gutorm/git/mobaGenetics-qc/qc-pipeline/snakefiles")
# print (f"DDDDDDDDEEEEBUG!  {snake_dir}   xxx")
# Grab the same configfiles as snakefile has
try: 
    with open(snake_dir/"rules.yaml", 'r') as stream:
        rule_info = yaml.safe_load(stream)
    with open(snake_dir/"config.yaml", 'r') as stream:
        config = yaml.safe_load(stream)
except Exception as e:
    print(f"mobaQcTools.py: Could load configfile, {str(e)}")

plink=config["plinklocal"]

def plotHist(dataFile, resultFile, column="name of the column", title="no legend??", separator='\s+',
             treshold=0, logx=False):
    """ plots and saves a historgram - can probably be removed as plot_point_and_line is now used

    Very basic Histogram. Could be prettied up a lot
    Prints out lots of warning, but see https://stackoverflow.com/questions/55805431/is-there-a-way-to-prevent-plotnine-from-printing-user-warnings-when-saving-ggplo
    Default separator is whitespace, but it needs to be overriden every now and then ... (typically by pure tab '\t')
    If (optional) logx is True, x-values are log10 transformed.
    """
    my_name = inspect.stack()[0][3]      # Generic way of find this functions name
    try: 
        df = pd.read_csv(dataFile, sep=separator, 
                 usecols=[column] )
    except Exception as e:
        print(f"{my_name}: Could not read plotdata {column} from {dataFile}, {str(e)}")
        return

    df = df.sort_values(column)
    p = ggplot(data=df, mapping=aes(x=column))
    hist = p + geom_histogram(binwidth=0.01) + labs(title=title)
    if treshold != 0:
        hist += geom_vline(xintercept=treshold, color='red')
    if logx:
        hist += scale_x_log10(name=f"log10({column})")
    ggsave(plot=hist, filename=resultFile, dpi=300)
    return

def plot_point_and_line(qc_results, dataFile, resultFile, column="name of the column", separator='\s+',
                     ylabel="no label given", invert=True):
    """ plots and saves a plot of sample/markers probabilities - probsabilities on y-axis

    Prints out lots of warning, but see https://stackoverflow.com/questions/55805431/is-there-a-way-to-prevent-plotnine-from-printing-user-warnings-when-saving-ggplo
    Default separator is whitespace, but it needs to be overriden every now and then ... (typically by pure tab '\t')
    Probabiliits are found as specified by dataFuke and column, with the separator given
    qc_results is a dictionary containing values that will be used for labels etc.
    If invert=True, will plot 1-probabilities . Default is True as this is often used to plot missingess/call rates.

    """
    my_name = inspect.stack()[0][3]      # Generic way of find this functions name. Use 
    try: 
        df = pd.read_csv(dataFile, sep=separator, 
                 usecols=[column] )
        treshold = qc_results.get("Treshold",0)
    except Exception as e:
        print(f"{my_name}: Could not read plotdata {column} from {dataFile} or qc-results, {str(e)}")
        return


    if invert:
        df = 1-df
        treshold = 1-treshold
    
    xlabel = qc_results.get("rule type")
    title = f'{qc_results.get("QC test")}\n{qc_results.get("Timestamp")}'
    df = df.sort_values(column).reset_index(drop=True)
    p = ggplot(df, aes(x=df.index,y=column)) 
    line = p + geom_line() + geom_point()
    if treshold > 0 and treshold <1 :
        line += geom_hline(yintercept=treshold, color='red')
        ylabel += f"   (treshold {treshold})"
        xlabel += f' ({qc_results.get("actionTakenCount")} outside treshold)' 
    line += labs(title=title, y=ylabel, x=xlabel)
    ggsave(plot=line, filename=resultFile, dpi=300)
    return

def saveYamlResults(files, yamlStruct):
    """ Creates three files from the given yamlStruct

    files is a Path or string, typically result.yaml
    "files.removedSamples" has yamlStruct["xitems"] intact
    "files" is the name of the structure without samples/markers (xitems). 
    A .rst file will be produced in addition to the .yaml file (used by snakemake --report)
    Side effect: Deletes yamlStruct["xitems"]
                 Sets Timestamp in yamlStruct
    """
    yamlStruct["Timestamp"] = datetime.utcnow().strftime("%Y.%m.%d %H:%M")+"(UTC)"
    files = Path(files)
    fullFile = files.with_suffix(".yaml.details")
    with open(fullFile, 'w') as file:
        yaml.dump(yamlStruct, file)
    # A shorter result version, without sample. 
    if "xitems" in yamlStruct:
        del yamlStruct["xitems"]
    with open(files, 'w') as file:
        yaml.dump(yamlStruct, file)
    # A .rst version used for captions
    rstFile = files.with_suffix(".rst")
    percentDropped = yamlStruct["actionTakenCount"] / yamlStruct["in"]
    with open(rstFile, 'w') as file:
        file.write(f'Rule {yamlStruct["Rule order"]} ({yamlStruct["rule type"]})\n\n')
        file.write(f'- {yamlStruct["in"]} in\n')
        file.write(f'- {yamlStruct["actionTakenCount"]} ({percentDropped:.1%}) {yamlStruct["rule action"]}d\n') # 
        file.write(f'- {yamlStruct["out"]} left\n\n')
        file.write(f'{yamlStruct["Timestamp"]}\n')

    
def plinkBase(path):
    """ part of the file without the last extention (such as .fam .bed .bim)

    plink often works on a trunk and creates extra files: a becomes a.fam, a.bed, a.bid
    This function reates this trunk, so /foo/bar/gazonk.fam becomes /foo/bar/gazonk
    The resulting trunk is typically used as input to plink
    Only the last '.' is removed so /foo.foo/bar/gazonk.x.fam will become /foo.foo/bar/gazonk.x
    """
    return  re.sub(r"\.\w*$","",path)


def lt(a,b):
    return a < b
def eq(a,b):
    return a == b
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


def extractSampleList(innFile, sampleFile, threshold_doc_file="/dev/null", 
                      colName="none", sep='\t', condition="<", treshold=0, cols={1}, ):
    """ A typical preprocessor to utilities like plink.

    Efficient in the sense that it 
    is not reading the (huge) files into memory.
    Takes a csv file innFile (first line with headers, seperator is paramneter "sep") and produces 
    * a sample list on sampleFile (one columun, for now hardcoded to column 0)
    * a threshold_doc_file with at set of numbered columns (cols) that also contains colName (see below)
    
    Only columns were the column colName matches the condition of treshold will be written
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
        subset = open(threshold_doc_file,"w+")
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

    

    
def dict_count_items(fil, cols=[0,1], warn=True):
    """
    Creates a dictionary with as key concatenation of the strings of the columns found in fil.
    Values are the number of times that key is found
    Returns the dictionary and the largest key value. In many scenarios, 1 is wanted here,
    meaning no duplicates. 
    if warn=True, will warn about key value > 1
    fil is expected to be a csv file with whitespace as delimiters
    First Column in the file is number 0 (standar Python)
    (Hint for when this is used to compare the columns of two files: You can reduce memory usage by 
    making a dictionary of the smallest file and searchin/iterating for matches through the largest

    Will print error/warning if dictionary is empty or if replicates found
    """
    # pandas is overkill here, but since we will use it anyway ...
    # Grab the relevant columns only
    try: 
        df = pd.read_csv(fil, usecols=cols, delim_whitespace=True, header=None)
    except Exception as e:
        print(f"Could not open file {fil}, {str(e)}")
        return
    # concat them as strings, handling NA
    df = pd.Series(df.fillna('NA').values.tolist()).map(lambda x: ' '.join(map(str,x)))
    # ... and finally make a dictionary. Note that we will here also count identical lines.
    countDict = dict()
    maxHits = 0
    for i in df:  
        countDict[i] = countDict.get(i,0) + 1
        if ( countDict[i]) > maxHits : 
            maxHits = countDict[i]
            sample = i
    if len(countDict) == 0 : print("ERROR: 0 sized dictionary after reading ",fil)
    if maxHits>1 and warn: print(f"WARNING: sample/marker {sample} found {maxHits} times in {fil}. Might not be the only nonunique ...")
    return countDict, maxHits

def lookupDict(fil, indx=1):
    """ Creates a lookup dictionary from fil, width column indx as index

    fil is a whitespice separated csv
    indx is the column that you will lookup on later, 0 is the first column
    """
    try: 
        all = pd.read_csv(fil, delim_whitespace=True, header=None).astype(str)
    except Exception as e:
        print(f"Could not open mapping file {fil}, {str(e)}")
        return
    indexCol = all[indx]  # get the index column
    # all = all.drop([indx], axis=1) # drop it - figured out it gave confusing output
    all = all.apply(" ".join, axis=1) # make each row a single string
    return dict(zip(indexCol,all))


def checkUpdates(preQc, postQc, cols=[0,1], indx=1, sanityCheck="none",
                 fullList=False,  mapFile="", mapIndx=1):
    """
    Return number of updates due to a QC-test as well as a structure suited for a export as a yaml-file
    The scenario is that a QC method/step had a dataset (file indData) and produced outData.
    Some items got filtered out/changed
    pre/postQC are tab-serparated csv-files with the same amount of columns
    The optional sanityCheck parameter will give error message as follows, depending to its value
       removal: Number of elements taken action on (removed) = size of preQc-postQc
       updated: size of preQc = postQc  
       anything different from the above): No tests performed
    Only columns passed by cols are used to compare input/and output.
    indx is only necessary if fullList is True. indx is used to genereate a list (usually of samples/markers)
    from column nr indx
    indx is sample-id/marker-id. indx=0 is the first column, default is second column. 
    mapFile is relevant when fullList=True and is the file that mapped values from preQc to postQc 
    and is used to show what the missing/renamed elements were named before
    mapIndx indicates the column of the mapFile that corresponds to indx in the datafile
    
    The yaml-structure will always contain the number of input samples/markers as well as sample/markers removed/remaining.
    
    """
    # dictionaly with only relevant columns
    (outDict,m) = dict_count_items(postQc, cols)
    result = {
        "in":   0,        # will count lines from the 'in' file
        "out": len(outDict), 
        "xitems": [],    # poplulated later by samples/markers if fullList is True
        "actionTakenCount": 0 # Not found in dictionary, so it is the effect of the qc-rule on the inputfile
        }
    haveDict = False
    if mapFile != "":     # Setting up a dictionary 'lookup' to 
                          # lookup the original value if postQc has changed
        lookup = lookupDict(mapFile, mapIndx)
        haveDict = True
        result["mapFile"] = mapFile
        #import json
        #print ("Saving lookupdict to lookup.json")
        #json.dump(lookup, open("lookup.json", 'w' ) )
        
    for line in open(preQc):
        result["in"] += 1 
        allcols = line.split()
        subsetc = [allcols[index] for index in cols]
        # concatenating so we can look up the strings with corresponding field from postQc file
        key = " ".join(map(str,subsetc))
        if (outDict.get(key,0) == 0):
            result["actionTakenCount"] += 1 
            item = str(allcols[indx])    # being the sample or the marker
            if fullList:
                changed = f"{item} - had {key}"  # This item was not found in preQc
                if haveDict:
                    mappingRule = lookup.get(item, 'ooops: missing')
                    changed = f'{changed} changed by mapping {mappingRule}'
                result["xitems"].append(changed)    # sample/marker id(s)
                
    if sanityCheck == 'updated':  # for updates, we dont't want to loose samples
        if (result["out"]) != result["in"] : 
            print(f"Warning: {preQc} -> {postQc}: Update expected, but number of unique items have changed "
                  f"from {result['in']} to {result['out']}" ) 
    elif sanityCheck == 'removal':       # for removals, we want out + removed = in
        if (result["actionTakenCount"] + result["out"]) != result["in"] : 
            print(f"Warning: {preQc} -> {postQc}: remaining + removed samples != original number")
    return result
    
def plinkBase(path):
    """ part of the file without the last extention (such as .fam .bed .bim)

    plink often works on a trunk and creates extra files: a becomes a.fam, a.bed, a.bid
    This function reates this trunk, so /foo/bar/gazonk.fam becomes /foo/bar/gazonk
    The resulting trunk is typically used as input to plink
    Only the last '.' is removed so /foo.foo/bar/gazonk.x.fam will become /foo.foo/bar/gazonk.x
    """
    return  re.sub(r"\.\w*$","",path)


def lt(a,b):
    return a < b
def eq(a,b):
    return a == b
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


def extractSampleList(innFile, sampleFile, threshold_doc_file = '/dev/null',
                      colName="none", condition="<", treshold=0, cols={1}, sep=None ):
    """ A typical preprocessor to utilities like plink.

    Despite the name, this also handles markers. Efficient in the sense that it 
    is not reading the (huge) files into memory.
    Takes a csv file innFile (first line with headers, seperator whitespace) and produces 
    * a sample list on sampleFile showing the numbered columns (cols)
    * a file whith the  numbered columns (cols) that also contains colName (see below)

    colName is a numerical column that will be compared to given threshold according to the
    wished threshold (using the condition parameter)
    For sample retrievel by plink, you will typically need the two first columns (fam/iid)
    while for marker extract you will typically need the first column (markerid)

    A seprator can be forced through sep, if not set one or more whitespaces is assumed to sperate columns
    
    Returns (number of sample extracted, total samples)
    Restrictions: 
    * Assumes only one match for the regExp colName
    * Columns are ordered in the same way as the innFile : {0,1} is identical to {1,0,0,1}
    * Outputfiles always have space as separator. They have no headers. 
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
        extra_info = open(threshold_doc_file,"w+")
        for line in enumerate(fp):
            allcols = line[1].split(sep)
            if compare(float(allcols[indx[0]]) , treshold): 
                # File without threshold
                subsetc = [allcols[index] for index in cols]
                sample.write(' '.join(map(str,subsetc)) + "\n" )
                # File with more columns, and always include the treshold column
                # subsetc is a relevant subset of the columns 
                extra_info.write(' '.join(map(str,subsetc)) + " " + allcols[indx[0]] + "\n" )
                matches += 1
    totalLines = line[0] + 1 # enumerates from 0, and we read a line manually
    sample.close()
    extra_info.close()
    return (matches,totalLines)

    

    
def dict_count_items(fil, cols=[0,1], warn=True):
    """
    Creates a dictionary with as key concatenation of the strings of the columns found in fil.
    Values are the number of times that key is found
    Returns the dictionary and the largest key value. In many scenarios, 1 is wanted here,
    meaning no duplicates. 
    if warn=True, will warn about key value > 1
    fil is expected to be a csv file with whitespace as delimiters
    First Column in the file is number 0 (standar Python)
    (Hint for when this is used to compare the columns of two files: You can reduce memory usage by 
    making a dictionary of the smallest file and searchin/iterating for matches through the largest

    Will print error/warning if dictionary is empty or if replicates found
    """
    # pandas is overkill here, but since we will use it anyway ...
    # Grab the relevant columns only
    try: 
        df = pd.read_csv(fil, usecols=cols, delim_whitespace=True, header=None)
    except Exception as e:
        print(f"Could not open file {fil}, {str(e)}")
        return
    # concat them as strings, handling NA
    df = pd.Series(df.fillna('NA').values.tolist()).map(lambda x: ' '.join(map(str,x)))
    # ... and finally make a dictionary. Note that we will here also count identical lines.
    countDict = dict()
    maxHits = 0
    for i in df:  
        countDict[i] = countDict.get(i,0) + 1
        if ( countDict[i]) > maxHits : 
            maxHits = countDict[i]
            sample = i
    if len(countDict) == 0 : print("ERROR: 0 sized dictionary after reading ",fil)
    if maxHits>1 and warn: print(f"WARNING: marker (or sample) identified by '{sample}' found {maxHits} times in {fil}. Might not be the only nonunique ...")
    return countDict, maxHits

def lookupDict(fil, indx=1):
    """ Creates a lookup dictionary from fil, width column indx as index

    fil is a whitespice separated csv
    indx is the column that you will lookup on later, 0 is the first column
    """
    try: 
        all = pd.read_csv(fil, delim_whitespace=True, header=None).astype(str)
    except Exception as e:
        print(f"Could not open mapping file {fil}, {str(e)}")
        return
    indexCol = all[indx]  # get the index column
    # all = all.drop([indx], axis=1) # drop it - figured out it gave confusing output
    all = all.apply(" ".join, axis=1) # make each row a single string
    return dict(zip(indexCol,all))


def checkUpdates(preQc, postQc, cols=[0,1], indx=1, sanityCheck="none",
                 fullList=False,  mapFile="", mapIndx=1):
    """
    Return number of updates due to a QC-test as well as a structure suited for a export as a yaml-file
    The scenario is that a QC method/step had a dataset (file indData) and produced outData.
    Some items got filtered out/changed
    pre/postQC are tab-serparated csv-files with the same amount of columns
    The optional sanityCheck parameter will give error message as follows, depending to its value
       removal: Number of elements taken action on (removed) = size of preQc-postQc
       updated: size of preQc = postQc  
       anything different from the above): No tests performed
    Only columns passed by cols are used to compare input/and output.
    fullList documents (in list xitems), what sample/markers have been changed/missing.
    indx is only necessary if fullList is True. indx is used to genereate a list (usually of samples/markers)
    from column nr indx
    indx is sample-id/marker-id. indx=0 is the first column, default is second column. 
    mapFile is relevant when fullList=True and is the file that mapped values from preQc to postQc 
    and is used to show what the missing/renamed elements were named before
    mapIndx indicates the column of the mapFile that corresponds to indx in the datafile
    
    The yaml-structure will always contain the number of input samples/markers as well as sample/markers removed/remaining.
    
    """
    # dictionaly with only relevant columns
    (outDict,m) = dict_count_items(postQc, cols)
    result = {
        "in":   0,        # will count lines from the 'in' file
        "out": len(outDict), 
        "xitems": [],    # poplulated later by samples/markers if fullList is True
        "actionTakenCount": 0 # Not found in dictionary, so it is the effect of the qc-rule on the inputfile
        }
    haveDict = False
    if mapFile != "":     # Setting up a dictionary 'lookup' to 
                          # lookup the original value if postQc has changed
        lookup = lookupDict(mapFile, mapIndx)
        haveDict = True
        result["mapFile"] = mapFile
        #import json
        #print ("Saving lookupdict to lookup.json")
        #json.dump(lookup, open("lookup.json", 'w' ) )
        
    for line in open(preQc):
        result["in"] += 1 
        allcols = line.split()
        subsetc = [allcols[index] for index in cols]
        # concatenating so we can look up the strings with corresponding field from postQc file
        key = " ".join(map(str,subsetc))
        if (outDict.get(key,0) == 0):
            result["actionTakenCount"] += 1 
            item = str(allcols[indx])    # being the sample or the marker
            if fullList:
                changed = f"{item} - had {key}"  # This item was not found in preQc
                if haveDict:
                    mappingRule = lookup.get(item, 'ooops: missing')
                    changed = f'{changed} changed by mapping {mappingRule}'
                result["xitems"].append(changed)    # sample/marker id(s)
                
    if sanityCheck == 'updated':  # for updates, we dont't want to loose samples
        if (result["out"]) != result["in"] : 
            print(f"Warning: {preQc} -> {postQc}: Update expected, but number of unique items have changed "
                  f"from {result['in']} to {result['out']}" ) 
    elif sanityCheck == 'removal':       # for removals, we want out + removed = in
        if (result["actionTakenCount"] + result["out"]) != result["in"] : 
            print(f"Warning: {preQc} -> {postQc}: remaining + removed samples != original number")
    return result


def log(logfile, message="Nothing logged", mode="a"):
    """
    Placeholder until we figure out how logging really needs to be done. Default appends message to logfile
    """
    with open(logfile, mode) as myfile:
        myfile.write(f'{datetime.utcnow().strftime("%Y%m%d %H:%M")} {message}')
    return

def _make_gen(reader):
    """
    Used to linecount buffered and fast
    """ 
    b = reader(1024*1024)
    while b:
        yield b
        b = reader(1024*1024)


def line_count(filename):
    """ Returns number of line if filename (assuming it is a text-file) just as wc -l

    Efficient linecounting by buffered reads
    """
    f = open(filename, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum( buf.count(b'\n') for buf in f_gen )


def create_exclude_list(duplicates, callRates, resultfile, excludelist):
    """
    Based on a file duplicates (created by plink --list-duplicate-vars), 
    and a call rate file (created by plink --missing)
    will create a list (resultfile) with the call rates/alleleles/etc of markers to be  excluded.
    Only the best call rate is kept, so the list will consist of the other duplicates
    Will also provide a simpler file (excludelist) suitable for plink to remove duplicates.
    If there were no duplicates, empty files are created.
    """
    if line_count(duplicates) == 1 :
        # Special case, no duplicates - only header in file duplicate
        mt = Path(resultfile)
        # From python 3.8, use unlink(missing_ok=True) . for now, live with the uglyness
        mt.touch()
        mt.unlink()
        mt.touch()
        mt = Path(excludelist)
        mt.touch()
        mt.unlink()
        mt.touch()
        return
    try: 
        calls = pd.read_csv(callRates, usecols=['CHR','SNP','F_MISS'],
                            delim_whitespace=True)
        dups = pd.read_csv(duplicates, usecols=['CHR','POS','ALLELES','IDS'], delimiter='\t')
    except Exception as e:
        print(f"Could not open inputfile, {str(e)}")        
        return

    # Series of ids, keep index - in order to explode the ids and join them back on chr
    ids = dups.IDS.str.split(' ',expand=True).stack().str.strip().reset_index(level=1, drop=True)
    df = pd.concat([dups.CHR, dups.POS, dups.ALLELES, ids], axis=1,
                   keys=['chr','pos','alleles','id'])
    df['duplicate_group'] = df.index   # The index show grouping, rememeber it
    # include call value to duplicate list
    all = df.merge(calls, left_on=['chr','id'], right_on=['CHR','SNP'])
                   
    # sort, group and keep all but the largest callrate
    drop=all.sort_values(by=['duplicate_group','F_MISS'],
             ascending=False).groupby('duplicate_group').apply(lambda g: g.iloc[1:])

    # Result file
    drop.loc[:,['chr', 'pos', 'alleles', 'SNP','F_MISS']].to_csv(resultfile,sep=" ", index=False)
    # snp idents plink laters can use
    drop.loc[:,['SNP']].to_csv(excludelist,sep=" ", index=False, header=False)
    print ("**** .... and we should make a yaml file with summary in it")
    

def exclude_strand_ambigious_markers(input, output, plink):
    """ Runs plink to exlucde A/T and C/G SNPs


    input is the trunk of a .bed/.bim/.fam triplet
    ouput will be the corresponding set, without excluded markers 
    The list over exluded markers can be found in ouput.excl
    plink will produce the output-bedset

    """
    try: 
        df = pd.read_csv(input+".bim", delim_whitespace=True, header=None).astype(str)
    except Exception as e:
        print(f"Could not open .bim-file of bedset {input}, {str(e)}")
        return
    mask = ((df[4] == 'G') & (df[5] == 'C') |
           (df[4] == 'C') & (df[5] == 'G') |
           (df[4] == 'A') & (df[5] == 'T') |
           (df[4] == 'T') & (df[5] == 'A'))
    (df[mask])[1].to_csv(output+".excl",index=False, header=False)

    subprocess.run([plink,
                "--bfile",input,
                "--exclude", output+".excl",                     
                "--out", output,               
                "--make-bed"
        ])


def exclude_strand_ambigious_markers(input, output, plink):
    """ Runs plink to exlucde A/T and C/G SNPs


    input is the trunk of a .bed/.bim/.fam triplet
    ouput will be the corresponding set, without excluded markers 
    The list over exluded markers can be found in ouput.excl
    plink will produce the output-bedset

    """
    try: 
        df = pd.read_csv(input+".bim", delim_whitespace=True, header=None).astype(str)
    except Exception as e:
        print(f"Could not open .bim-file of bedset {input}, {str(e)}")
        return
    mask = ((df[4] == 'G') & (df[5] == 'C') |
           (df[4] == 'C') & (df[5] == 'G') |
           (df[4] == 'A') & (df[5] == 'T') |
           (df[4] == 'T') & (df[5] == 'A'))
    (df[mask])[1].to_csv(output+".excl",index=False, header=False)

    subprocess.run([plink,
                "--bfile",input,
                "--exclude", output+".excl",                     
                "--out", output,               
                "--make-bed"
        ])
    
    
def missing_genotype_rate(rule,
                          in_bedset, out_bedset, sample=True, treshold=0.1,
                          result_file='/dev/null', plot_file=False):
    """ Runs plink --geno or --mind and produces output

    Wrapper around plink and saveYamlResults as well as plots showing what gets removed
    If plot_file is False, no plot is produced
    Returns a 'dropouts' structure that contains info about the dropped samples/markers 
    (but only summaries).
    rule is the calling rule, and will be added to the result_file/dropouts
    
    The bedset are truncs and not files: For example for in_bedset=foo and sample=True, 
    plink --mind will be used, and exclusion results will be related to foo.fam
    
    """
    if sample:
        kindof = "sample"
        extension = ".fam"
        plink_switch = "--mind"
        miss_ext = ".imiss"   # if plotting, .imiss is for samples, .lmiss for markers
    else: # marker
        kindof = "marker"
        extension = ".bim"
        plink_switch = "--geno"
        miss_ext = ".lmiss"
    
    subprocess.run([plink,
                "--bfile",in_bedset,
                plink_switch, str(treshold),
                "--out", out_bedset,               
                "--make-bed"
        ])
    
    dropouts = checkUpdates(in_bedset+extension, out_bedset+extension, cols = [0,1],
                            sanityCheck = "removal", fullList = True) 
    dropouts.update(rule_info[rule])   # Metainfo and documentation about the rule
    dropouts["Treshold"] = treshold
    dropouts["Rule"] = rule
    dropouts["rule type"] = kindof # rule says sample/marker - we know what is really is
    saveYamlResults(result_file, dropouts)
    # A plot might be ordered
    if plot_file:
        # call rates for markers/samples before they got removed
        subprocess.run([plink,
                "--bfile", in_bedset,
                "--missing",
                "--out", in_bedset ])
        plot_point_and_line(dropouts, in_bedset+miss_ext, plot_file,
                                column="F_MISS",ylabel="1 - missingness")

    return dropouts

def low_hwe_autosomal_rate(rule,
                          in_bedset, out_bedset, sample=True, treshold=0.1, #remove sample?
                          result_file='/dev/null', plot_file=False):
    """ Runs plink plink autosomal produces output, including plot

    Wrapper around plink and saveYamlResults as well as plots showing what gets removed
    If plot_file is False, no plot is produced
    Returns a 'dropouts' structure that contains info about the dropped samples/markers 
    (but only summaries).
    rule is the calling rule, and will be added to the result_file/dropouts
    
    The bedset are truncs and not files: For example for in_bedset=foo and sample=True, 
    plink --mind will be used, and exclusion results will be related to foo.fam
    
    """
    print ("DEBUGGGGGGGGGGGGG Clean up sample mess")
    if sample:
        kindof = "sample"
        extension = ".fam"
        plink_switch = "--mind"
        miss_ext = ".imiss"   # if plotting, .imiss is for samples, .lmiss for markers
    else: # marker
        kindof = "marker"
        extension = ".bim"
        plink_switch = "--geno"
        miss_ext = ".lmiss"
    
    subprocess.run([plink,
                "--bfile",in_bedset,
                "--autosome",
                "--hardy", "midp",
                "--out", out_bedset,               
        ])
    # We here have a .hwe file where low p-falues are to be removed
    extractSampleList(out_bedset+".hwe", out_bedset+".exclude",
                      threshold_doc_file=out_bedset+".details",
                      colName="^P$", condition="<", treshold=treshold, cols={1})

    
    # dropouts = checkUpdates(in_bedset+extension, out_bedset+extension, cols = [0,1],
    #                         sanityCheck = "removal", fullList = True) 
    # dropouts.update(rule_info[rule])   # Metainfo and documentation about the rule
    # dropouts["Treshold"] = treshold
    # dropouts["Rule"] = rule
    # dropouts["rule type"] = kindof # rule says sample/marker - we know what is really is
    # saveYamlResults(result_file, dropouts)
    # # A plot might be ordered
    # if plot_file:
    #     # call rates for markers/samples before they got removed
    #     subprocess.run([plink,
    #             "--bfile", in_bedset,
    #             "--missing",
    #             "--out", in_bedset ])
    #     plot_point_and_line(dropouts, in_bedset+miss_ext, plot_file,
    #                             column="F_MISS",ylabel="1 - missingness")

#    return dropouts

def exclude_strand_ambigious_markers(input, output, plink):
    """ Runs plink to exlucde A/T and C/G SNPs
    input is the trunk of a .bed/.bim/.fam triplet
    ouput will be the corresponding set, without excluded markers 
    The list over exluded markers can be found in ouput.excl
    plink will produce the output-bedset
    """
    try: 
        df = pd.read_csv(input+".bim", delim_whitespace=True, header=None).astype(str)
    except Exception as e:
        print(f"Could not open .bim-file of bedset {input}, {str(e)}")
        return
    mask = ((df[4] == 'G') & (df[5] == 'C') |
           (df[4] == 'C') & (df[5] == 'G') |
           (df[4] == 'A') & (df[5] == 'T') |
           (df[4] == 'T') & (df[5] == 'A'))
    (df[mask])[1].to_csv(output+".excl",index=False, header=False)

    subprocess.run([plink,
                "--bfile",input,
                "--exclude", output+".excl",                     
                "--out", output,               
                "--make-bed"
        ])
def fix_rsid_map(mapfile, newmap):
    """ Create a rsid mapping based on som Moba business-logic
    Map of the rsid given in mapfile needs tweeking. newmap is produced
    Multiple whitespaces in mapfile will become a single space in newmap
    Could have been a lot more efficient
    We do the following: (works for GSA/GSADM)
    * Ignore lines that map to dot ('.')
    * If multiple (comma-separated) ids are found in to-map, we use the first
    * from strings containg .1 .2 ... .9 are ignored
    """
    try: 
        mappings = pd.read_csv(mapfile, usecols=[0,1], names=['from','to'], 
                            delim_whitespace=True).astype(str)
    except Exception as e:
        print(f"Could not open file {mapfile}, {str(e)}")        
        return

    with open(newmap, "w") as out:
        multiAllele = re.compile("\.\d")
        for index,row in mappings.iterrows():
            # Elements to ignore
            if row['to'] == "." : continue
            if multiAllele.search(row['from']) : continue   #eg rs222.1 is ignored
            # Elements to simplify:
            # 1 - Trucate everything including and after the first comma
            to = re.sub(r",.+$","",row['to'])
            # Save
            out.write(f"{row['from']} {to}\n")


def intersect_rsid(bim_small, bim_big, intersection):
    """ Assumes bim files, that is tab-serarated plink with rsID in second column

    intersection is a file to be created
    If one of the files are large, pass that as bim_big for efficiency
    
    could have been more general checking any column

    """
    my_name = inspect.stack()[0][3]      # Generic way of find this functions name
    m = 0   # max hits
    try: 
        (smallDict,m) = dict_count_items(bim_small, [1], warn=False) # what do we have
        with open(intersection, "w") as out:
            for line in open(bim_big):
                rsid = line.split()[1]
                if smallDict.get(rsid,0) > 0 :  out.write(f"{rsid}\n")

    except Exception as e:
        print(f"{my_name} Exception caught:  {str(e)}")        
        
            
def copy_file(f,ext=".bak"):
    """ makes a .bak file

    This is a debug-tool used keep a copy of a snakemake result that it would have deleted due 
    to failure.

    """
    print(f"DEBUG> Making a backup of {f} to {f+ext}.")
    copyfile(f, f+ext)

def dotplot(genomedata,prec=2,x='x',y='y',c='c'):
    """ Returns a plotnine object ready to be printet, but where extra lines can be added

    Assumes x and y are numbers, these will be rounded to precision decimals
    The preicision is there to not cluster set plot when there are two many points to be seen
    c is used for colouring and shaping - a colourblind palette will be used
    The data is found in a whitespace separated file where they x,y and c are headers
    
    """
    my_name = inspect.stack()[0][3]      # Generic way of find this functions name        
    try: 
        df = pd.read_csv(genomedata, usecols=[c,x,y], delim_whitespace=True)
    except Exception as e:
        print(f"{my_name}: {str(e)}")        
        return
    
    p = ggplot(data=df.round(prec).drop_duplicates(), mapping=aes(x=x,y=y,color=c, shape=c))
    p += scale_colour_brewer(type="qual", palette="Set1")  # better for colourblind
    p += geom_point()
    
    return p


