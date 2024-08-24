#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import math
import plotnine as p9
# While plotting in a non-interactive environment (xterm ...) we need this:
import yaml
import re
from pathlib import Path
import subprocess
from datetime import datetime
from shutil import copyfile
import inspect    # to find stack-info, such as the functions name
import matplotlib
import gzip
import shutil
matplotlib.use('Agg')


def make_rule_caption(rule_info, rule, dir):
    """ Creates a caption for rule

    Rule is the calling rule, dir is a path object where the file rule.rst is to be created

    Note that this is a hack because founder/offspring construction for now (feb 2020) has
    proven difficutlt to make good role-dependent caption for.
    Since captions are used for sorting the html report produced by snakemake --report,
    this function creates one based on the global rule_info[] structure.
    It will typically be called twice.
    Warning: Assumes that rule and rule_stem are identical where this is used.
    """
    with open(dir/(rule + ".rst"), 'w') as file:
        file.write(f'Rule {rule_info[rule]["Rule order"]} ({rule_info[rule]["rule action"]})\n')


def plot_hist(dataFile, resultFile, column="name of the column",
              title="no legend??", separator='\s+',
              threshold=0, logx=False, bins=100):
    """ plots and saves a histogram

    Prints out lots of warning, but see https://stackoverflow.com/questions/55805431/is-there-a-way-to-prevent-plotnine-from-printing-user-warnings-when-saving-ggplo
    Default separator is whitespace, but it needs to be overriden every now and then ... (typically by pure tab '\t')
    If (optional) logx is True, x-values are log10 transformed.
    """
    my_name = inspect.currentframe().f_code.co_name  # Generic way of find function name
    try:
        df = pd.read_csv(dataFile, sep=separator, usecols=[column])
    except Exception as e:
        print(f"{my_name}: Could not read plotdata {column} from {dataFile}, {str(e)}")
        return

    df = df.sort_values(column)
    p = p9.ggplot(
        data = df,
        mapping = p9.aes(
            x = column
        )
    )

    hist = p + p9.geom_histogram(
        bins = bins
    )
    hist += p9.labs(
        title = title,
        x = column
    )
    if threshold != 0:
        hist += p9.geom_vline(
            xintercept = threshold,
            color = 'red'
        )
    if logx:
        hist += p9.scale_x_log10(
            name = f"log10({column})"
        )
    hist.save(
        filename = resultFile,
        dpi = 300
    )
    return


def plot_point_and_line(
        qc_results,
        dataFile,
        resultFile,
        column = "name of the column",
        separator = '\s+',
        ylabel = "no label given",
        invert = True
):
    """plots and saves a plot of sample/markers probabilities (y-axis)

    Prints out lots of warning, but see
    https://stackoverflow.com/questions/55805431/is-there-a-way-to-prevent-plotnine-from-printing-user-warnings-when-saving-ggplo
    Default separator is whitespace, but it needs to be overriden
    every now and then ... (typically by pure tab '\t') Probabiliits
    are found as specified by dataFile and column, with the separator
    given qc_results is a dictionary containing values that will be
    used for labels etc.  If invert=True, will plot 1-probabilities.
    Default is True as this is often used to plot missingess/call
    rates.

    """
    my_name = inspect.currentframe().f_code.co_name
    try:
        df = pd.read_csv(
            dataFile,
            sep = separator,
            usecols = [column]
        )
        threshold = qc_results.get("Threshold", 0)
    except Exception as e:
        print(f"{my_name}: Could not read plotdata {column} from {dataFile} or qc-results, {str(e)}")
        return

    if invert:
        df = 1-df
        threshold = 1-threshold

    xlabel = qc_results.get("rule type")
    title = f'{qc_results.get("QC test")}\nthreshold={threshold}'
    df = df.sort_values(column).reset_index(drop=True)

    p = p9.ggplot(
        df,
        p9.aes(
            x = df.index,
            y = column
        )
    )

    line = p + p9.geom_line() + p9.geom_point()

    if threshold > 0 and threshold < 1 :
        line += p9.geom_hline(
            yintercept = threshold,
            color = 'red'
        )
        # ylabel += f"   (threshold={threshold})"
        xlabel += f' ({qc_results.get("actionTakenCount")} outside threshold)'

    line += p9.labs(
        title = title,
        y = ylabel,
        x = xlabel
    )

    print("- Exporting to ", resultFile)

    line.save(
        filename = resultFile,
        dpi = 300,
        width = 4,
        height = 3,
        units = "cm"
    )
    return


def plot_text(
        text,
        plotFile
):
    d = {
        'x': [0],
        'y': [0],
        'label': [text]
    }

    df = pd.DataFrame(data=d)

    p = p9.ggplot() + p9.geom_text(
        df,
        p9.aes(
            x = 'x',
            y = 'y',
            label = 'label'
        )
    )

    p.save(
        filename = plotFile,
        dpi = 300,
        width = 4,
        height = 3,
        units = "cm"
    )
    return


def saveYamlResults(files, yamlStruct):
    """ Creates 3 files (.rst .yaml .yaml.details) from yamlStruct

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
    # March 2020; Note that founder/offspring rules, it has been hard
    # to declare separate rst- files for founder/offspring. Instead,
    # we later create a common file, overwriting whatever is made here
    # :-(
    with open(rstFile, 'w') as file:
        file.write(f'Rule {yamlStruct["Rule order"]} ({yamlStruct["rule type"]})\n\n')
        file.write(f'- {yamlStruct["in"]} in\n')
        file.write(f'- {yamlStruct["actionTakenCount"]} ({percentDropped:.1%}) {yamlStruct["rule action"]}\n') #
        file.write(f'- {yamlStruct["out"]} left\n\n')
        file.write(f'{yamlStruct["Timestamp"]}\n')


def plinkBase(path):
    """ part of the file without the last extension (such as .fam .bed .bim)

    plink often works on a trunk and creates extra files: a becomes a.fam, a.bed, a.bid
    This function creates this trunk, so /foo/bar/gazonk.fam becomes /foo/bar/gazonk
    The resulting trunk is typically used as input to plink
    Only the last '.' is removed so /foo.foo/bar/gazonk.x.fam will become /foo.foo/bar/gazonk.x
    """
    return re.sub(r"\.\w*$", "", path)


def lt(a, b):
    try:
        result = (a < b)
    except Exception as e:
        print(f"Could not compare '{a}' < '{b}', {str(e)}")
    return result


def eq(a, b):
    try:
        result = (a == b)
    except Exception as e:
        print(f"Could not compare '{a}' == '{b}', {str(e)}")
    return result


def gt(a, b):
    try:
        result = (a > b)
    except Exception as e:
        print(f"Could not compare '{a}' > '{b}', {str(e)}")
    return result


def unknownComp(a, b):
    print("Unknown compare operation")
    return


# used to map symbols to actual comparision functions defined previously
compOperands = {
    "<": lt,
    "==": eq,
    ">": gt
}


def extract_list(
        inFile,
        outFile,
        threshold_doc_file = "/dev/null",
        colName = "none",
        sep = '\t',
        condition = "<",
        threshold = 0,
        key_cols = [0],
        doc_cols = [0, 1]
):
    """A typical preprocessor to utilities like plink, extracts relevant samples/markers

    Efficient in the sense that it is not reading the (huge) files
    into memory.  Takes a csv file inFile (first line with headers,
    separator is parameter "sep", use None for whitespaces)

    Produces
    * a  list on outFile of the columns indicated by key_cols.
    * a threshold_doc_file with at set of numbered columns (doc_cols).
      The contents of colName will always be added

    Only columns were the column colName (regexp) matches the
    condition of threshold will be written threshold can be a string,
    this is testet only for condition '==' Returns (number of sample
    extracted, total samples)

    Restrictions:
    * Prints error and returns unless only one columns matches
    * Outputfiles use the same separator as used for inFile. If 'none'
     was used ' ' is used

    """
    print(inFile)
    with open(inFile) as fp:

        # Identifying header column
        line = fp.readline()
        allcols = line.split(sep)
        outsep = sep

        if outsep is None:
            outsep = " "

        regex = re.compile(colName)
        indx = [i for i, item in enumerate(allcols) if re.search(regex, item)]

        if len(indx) != 1:
            print(f"ERROR in {inspect.currentframe().f_code.co_name}: regexp {colName}  found {len(indx)} times in {inFile}. First line was {allcols}")
            return

        indx = indx[0]

        matches = 0
        compare = compOperands.get(condition, unknownComp)
        if compare == unknownComp:
            print(f"Cannot compare using: '{condition}'")
            return

        sample = open(outFile, "w+")
        subset = open(threshold_doc_file, "w+")

        for line in enumerate(fp):

            allcols = line[1].split(sep)
            val = allcols[indx]
            try: val = float(val)  # Skip non-numeric values
            except: pass

            if compare(val, threshold):

                # File with only keys, create the relevant subset first
                subsetc = [allcols[index] for index in key_cols]
                sample.write(f"{outsep.join(map(str, subsetc))}\n")
                # File with more columns, and always include the threshold column
                # subsetc is a relevant subset of the columns
                subsetc = [allcols[index] for index in doc_cols]
                subset.write(f"{outsep.join(map(str, subsetc))} {val}\n")
                matches += 1

    totalLines = line[0] + 1  # enumerates from 0, and we read a line manually
    sample.close()
    subset.close()

    return (matches, totalLines)


def dict_count_items(fil, cols=[0, 1], warn=True):
    """Creates and counts items in a dictionary

    Creates a dictionary with as key concatenation of the strings of
    the columns cols found in fil.
    Values are the number of times that key is found
    Returns the dictionary and the largest key value. In many
    scenarios, 1 is wanted here, meaning no duplicates.
    if warn=True, will warn about key value > 1
    fil is expected to be  csv file with whitespace as delimiters
    First Column in the file is number 0 (standard Python)
    (Hint for when this is used to compare the columns of two files:
    You can reduce memory usage by making a dictionary of the smallest
    file and searchin/iterating for matches through the largest

    Will print error/warning if dictionary is empty or if replicates found

    """
    # Grab the relevant columns only
    try:
        df = pd.read_csv(
            fil,
            usecols = cols,
            sep = '\s+',
            header = None
        )
    except Exception as e:
        print(f"Could not open file {fil}, {str(e)}")
        return

    # concat them as strings, handling NA
    df = pd.Series(df.fillna('NA').values.tolist()).map(lambda x: ' '.join(map(str, x)))
    # ... and finally make a dictionary. Note that we will here also count identical lines.
    countDict = dict()
    maxHits = 0
    for i in df:
        countDict[i] = countDict.get(i, 0) + 1
        if (countDict[i]) > maxHits:
            maxHits = countDict[i]
            sample = i
    if len(countDict) == 0:
        print("ERROR: 0 sized dictionary after reading ",fil)
    if maxHits > 1 and warn:
        print(f"WARNING: sample/marker {sample} found {maxHits} times in {fil}. Might not be the only nonunique ...")
    # print(f"**********   counted {fil} and got {len(countDict)} unique items")
    return countDict, maxHits


def lookupDict(fil, indx=1):
    """ Creates a lookup dictionary from fil, width column indx as index

    fil is a whitespice separated csv
    indx is the column that you will lookup on later, 0 is the first column
    """
    try:
        all = pd.read_csv(
            fil,
        sep = '\s+',
            header = None
        ).astype(str)
    except Exception as e:
        print(f"Could not open mapping file {fil}, {str(e)}")
        return
    indexCol = all[indx]  # get the index column
    # all = all.drop([indx], axis=1)   # drop it - figured out it gave confusing output
    all = all.apply(" ".join, axis=1)  # make each row a single string
    return dict(zip(indexCol, all))


def checkUpdates(
        preQc,
        postQc,
        cols = [0, 1],
        indx = 1,
        sanityCheck = "none",
        fullList = False,
        mapFile = "",
        mapIndx = 1,
        allele_flip = False
):
    """Generic comparition of two bedsets

    Return number of updates due (typically) to a rule as well as a
    structure suited for a export as a yaml-file

    The scenario is that a method/step had a dataset (file preQc)
    and produced postQc.

    Some items got filtered out/changed. pre/postQC are tab-serparated
    csv-files with the same amount of columns

    The optional sanityCheck parameter will give error message as
    follows, depending to its value
       removal: Number of action on (removed) = size of preQc-postQc
       updated: size of preQc = postQc
       anything different from the above): No tests performed
    Only columns passed by cols are used to compare between input/and output.

    fullList documents (in list xitems), what sample/markers have been
    changed/missing.  

    indx (if fullList is True) is used to genereate a list (usually of
    samples/markers) from column nr indx.
    indx is sample-id/marker-id. indx=0 is the first column, default
    is second column.

    mapFile is relevant when fullList=True and is the file that mapped
    values from preQc to postQc and is used to show what the
    missing/renamed elements were named before

    mapIndx indicates the column of the mapFile that corresponds to
    indx in the datafile

    if allele_flip is true, the two last columns are expected to be
    letters (alleles). Records will be classified as identical
    regardless of the order of these: "foo A T" will be regarded as
    identical too "foo T A" (but obiously not of "foo T x" etc)

    The yaml-structure will always contain the number of input
    samples/markers as well as sample/markers removed/remaining.

    """
    # dictionary with only relevant columns
    (outDict, m) = dict_count_items(postQc, cols)
    result = {
        "in":   0,        # will count lines from the 'in' file
        "out": len(outDict),
        "xitems": [],    # populated later by samples/markers if fullList is True
        "actionTakenCount": 0 # Not found in dictionary, so it is the effect of the qc-rule on the inputfile
    }
    haveDict = False
    if mapFile != "":   # Setting up a dictionary 'lookup' to look up the original value if postQc has changed
        lookup = lookupDict(mapFile, mapIndx)
        haveDict = True
        result["mapFile"] = mapFile
    for line in open(preQc):
        result["in"] += 1
        allcols = line.split()
        subsetc = [allcols[index] for index in cols]
        # concatenating so we can look up the strings with corresponding field from postQc file
        key = " ".join(map(str, subsetc))
        flipkey = key  # in case we want to flip alleles, if not key == flipkey
        if allele_flip:
            parts = re.match(r"(.*) (\S+) (\S+)$", key) # assuming the alleles at the end
            # group(2) is allele 1, group(3) allele 2
            # Turned out that alleles (like in 1000genomes can be on the form
            # <INS:ME:ALU> so we are accept any non-blank in the match
            #print (key)
            flipkey = " ".join((parts.group(1),parts.group(3),parts.group(2)))
        if (outDict.get(key, 0) + outDict.get(flipkey, 0)) == 0 :
            # If we didnt find this, even by (possible) flipping alleles
            result["actionTakenCount"] += 1
            item = str(allcols[indx])    # being the sample or the marker
            if fullList:
                changed = f"{item} - had {key}"  # This item was not found in preQc
                if haveDict:
                    mappingRule = lookup.get(item, 'ooops: missing')
                    changed = f'{changed} changed by mapping {mappingRule}'
                result["xitems"].append(changed)    # sample/marker id(s)

    if sanityCheck == 'updated':  # we don't want to lose items here
        if (result["out"]) != result["in"]:
            print(f"Warning: {preQc} -> {postQc}: Update expected, but number of unique items have changed "
                  f"from {result['in']} to {result['out']}")
    elif sanityCheck == 'removal':       # we want out + removed = in
        if (result["actionTakenCount"] + result["out"]) != result["in"]:
            print(f"Warning: {preQc} -> {postQc}: remaining + removed samples != original number")

    return result


def log(logfile, message="Nothing logged", mode="a"):
    """ Placeholder until we figure out how logging really needs to be
    done. Default appends message to logfile

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
    """ Returns number of lines is (text) filename as wc -l

    Efficient linecounting by buffered reads
    """
    f = open(filename, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum( buf.count(b'\n') for buf in f_gen)


def create_exclude_list(duplicates, callRates, resultfile, excludelist):
    """Based on a file duplicates (created by plink --list-duplicate-vars), and a call rate file (created by plink --missing)
    will create a list (resultfile) with the call rates/alleleles/etc
    of markers to be excluded.
    Only the best call rate is kept, so the list will consist of the
    other duplicates
    Will also provide a simpler file (excludelist) suitable for plink
    to remove duplicates.  If there were no duplicates, empty files
    are created.

    """
    if line_count(duplicates) == 1:
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
        calls = pd.read_csv(
            callRates,
            usecols = ['CHR', 'SNP', 'F_MISS'],
            sep = '\s+'
        )
        dups = pd.read_csv(
            duplicates,
            usecols = ['CHR', 'POS', 'ALLELES', 'IDS'],
            delimiter='\t'
        )
    except Exception as e:
        print(f"Could not open inputfile, {str(e)}")
        return

    # Series of ids, keep index
    # in order to explode the ids and join them back on chr
    ids = dups.IDS.str.split(' ',expand=True).stack().str.strip().reset_index(level=1, drop=True)
    df = pd.concat([dups.CHR, dups.POS, dups.ALLELES, ids], axis=1,
                   keys=['chr', 'pos', 'alleles', 'id'])
    df['duplicate_group'] = df.index   # The index show grouping, rememeber it
    # include call value to duplicate list
    all = df.merge(calls, left_on = ['chr', 'id'], right_on = ['CHR', 'SNP'])

    # sort, group and keep all but the largest callrate
    drop = all.sort_values(
        by = ['duplicate_group', 'F_MISS'],
        ascending = False
    ).groupby(
        'duplicate_group'
    ).apply(
        lambda g: g.iloc[1:]
    )

    # Result file
    drop.loc[:, ['chr', 'pos', 'alleles', 'SNP','F_MISS']].to_csv(resultfile, sep = " ", index = False)
    # snp idents plink laters can use
    drop.loc[:, ['SNP']].to_csv(excludelist, sep=" ", index = False, header = False)
    print("**** .... and we should make a yaml file with summary in it")


def missing_genotype_rate(
        rule,
        in_bedset,
        out_bedset,
        sample = True,
        threshold = 0.1,
        result_file = '/dev/null',
        plot_file = False,
        plinklocal = None,
        rule_info = None
):
    """Runs plink --geno or --mind and produces output

    Wrapper around plink and saveYamlResults as well as plots showing
    what gets removed

    If plot_file is False, no plot is produced

    Returns a 'dropouts' structure that contains info about the
    dropped samples/markers (but only summaries). rule is the calling
    rule, and will be added to the result_file/dropouts

    The bedset are truncs and not files: For example for in_bedset=foo
    and sample=True, plink --mind will be used, and exclusion results
    will be related to foo.fam

    """
    if sample:
        plink_switch = "--mind"
    else:  # marker
        plink_switch = "--geno"


    subprocess.run(
        [
            plinklocal,
            "--bfile", in_bedset,
            plink_switch, str(threshold),
            "--out", out_bedset,
            "--make-bed"
        ],
        check = True
    )

    return missing_genotype_rate_docs(
        rule = rule,
        in_bedset = in_bedset,
        out_bedset = out_bedset,
        sample = sample,
        threshold = threshold,
        result_file = result_file,
        plot_file = plot_file,
        plinklocal = plinklocal,
        rule_info = rule_info
    )

def missing_genotype_rate_docs(
        rule,
        in_bedset,
        out_bedset,
        sample = True,
        threshold = 0.1,
        result_file = '/dev/null',
        plot_file = False,
        plinklocal = None,
        rule_info = None
):
    """Computes dropout information and makes plots related to a missing_genotype_rate call.

    """

    if sample:
        kindof = "sample"
        extension = ".fam"
        miss_ext = ".imiss"   # if plotting, .imiss is for samples, .lmiss for markers
    else:  # marker
        kindof = "marker"
        extension = ".bim"
        miss_ext = ".lmiss"

    dropouts = checkUpdates(
        str(in_bedset) + extension,
        str(out_bedset) + extension,
        cols = [0, 1],
        sanityCheck = "removal",
        fullList = True
    )

    dropouts.update(rule_info[rule])   # Metainfo and documentation about the rule
    dropouts["Threshold"] = threshold
    dropouts["Rule"] = rule
    dropouts["rule type"] = kindof  # rule says sample/marker - we know what it really is
    saveYamlResults(result_file, dropouts)

    if plot_file:
        # call rates for markers/samples before they got removed
        subprocess.run(
            [
                plinklocal,
                "--bfile", in_bedset,
                "--missing",
                "--out", in_bedset
            ],
            check = True
        )
        plot_point_and_line(
            dropouts,
            str(in_bedset) + miss_ext,
            plot_file,
            column = "F_MISS",
            ylabel = "1 - missingness"
        )

    return dropouts


def compute_hwe(
        in_bedset,
        out_bedset,
        threshold = 0.1,
        hwe_switches = ["--autosome", "--hardy", "midp"],
        plinklocal = None
):
    """Runs plink hwe and computes p-values but does not change .bed file

    P-values for markers will end up in out_bedset.hwe - a file
    suitable for extract_list but not for plink.
    A subset below threshold will be found in out_bedset.exclude
    (suitable for plink exclude) and out_bedset.details will contain
    details including p-values

    """
    subprocess.run(
        [
            plinklocal,
            "--bfile", in_bedset,
            "--out", out_bedset
        ] + hwe_switches,
        check = True
    )

    # We here have a .hwe file where low p-values for markers are to be removed
    hwe_p_values = out_bedset + ".hwe"

    extract_list(
        hwe_p_values,
        out_bedset + ".exclude",
        threshold_doc_file = out_bedset + ".details",
        sep = None,
        colName = "^P$",
        condition = "<",
        threshold = threshold,
        key_cols = [1],
        doc_cols = [0, 1]
    )

    return

def filter_hwe(
        rule,
        in_bedset,
        out_bedset,
        threshold = 0.1,
        hwe_switches = ["--autosome", "--hardy", "midp"],
        result_file = '/dev/null',
        plot_file = '/dev/null',
        plinklocal = None,
        rule_info = None
):
    """ Runs plink hwe and removes low p-values markers. Produces output, including plot.

    Wrapper around plink to do Hardy Weinberg Equilibrium tests and plot distribution.
    Saves results with  saveYamlResults as well.
    bedsets are 'trunks'

    """

    compute_hwe(
        in_bedset,
        out_bedset,
        threshold,
        hwe_switches,
        plinklocal = plinklocal
    )

    hwe_p_values = out_bedset + ".hwe"

    # low values detected, now extract, make results and plot
    subprocess.run(
        [
            plinklocal,
            "--bfile", in_bedset,
            "--exclude", out_bedset + ".exclude",
            "--out", out_bedset,
            "--make-bed"
        ],
        check = True
    )

    dropouts = checkUpdates(
        in_bedset + ".bim",
        out_bedset + ".bim",
        cols = [0, 1],
        sanityCheck = "removal",
        fullList = True
    )

    dropouts.update(rule_info[rule])   # Metainfo and documentation about the rule
    dropouts["Threshold"] = threshold
    dropouts["Rule"] = rule
    saveYamlResults(result_file, dropouts)

    p = hweg_qq_plot(hwe_p_values, prec = 2, x = 'P')
    t_line = p9.geom_hline(
        yintercept = -1 * math.log10(threshold),
        color = 'red'
    )
    p += t_line    # threshold. Note that plot removed above
    title = f'HWE > -log({threshold}) \n{dropouts.get("actionTakenCount")} outside threshold\n{dropouts.get("Timestamp")}'
    p += p9.labs(title = title)
    p.save(
        filename = plot_file,
        dpi = 600
    )
    return


def compute_excess_het(het_file, out_file, sd):
    """Computes a list of samples with excess het

    Based on het_file and a number of standard deviation sd,
    Creates three files:
         out_file contains samples (family, iid) to be removed becase
         they have --het excess too large (comparing het_rate (se
         code) to sd*standard deviations) .

         This file is intented to be used by plink
         to remove these samples and has no headers

         out_files.details contains information on how they actually
         were removed. This file has headers

         out_files.total contains .het file + more metrics that caused
         remival. This file has headers

    """
    my_name = inspect.currentframe().f_code.co_name      # This functions name
    try:
        df = pd.read_csv(
            het_file,
            sep = '\s+'
        )
    except Exception as e:
        print(f"{my_name}: {str(e)}")
        return

    # Columns defined in https://www.cog-genomics.org/plink/1.9/formats#het
    df['het_rate'] = (df['N(NM)'] - df['O(HOM)']) / df['N(NM)']
    failed = df[df.het_rate > df.het_rate.mean() + sd*df.het_rate.std()].copy()
    failed['het_dst'] = (failed['het_rate'] - df['het_rate'].mean()) / df['het_rate'].std()

    # output - the .total file is used for plotting later
    df.to_csv(
        out_file + ".total",
        sep = " ",
        index = False
    )

    # file to be used by plink to remove samples
    failed[['FID', 'IID']].to_csv(
        out_file,
        sep = " ",
        index = False,
        header = False
    )

    # Documentation of the actual metrics  of the removed samples
    failed.to_csv(
        out_file + ".details",
        sep = " ",
        index = False
    )

def extract_gz(
        gzfile,
        destinationFile
):
    """
    Extract the content of a gzipped file to another file
    """

    with gzip.open(gzfile, 'rb') as file_in:
        with open(destinationFile, 'wb') as file_out:
            shutil.copyfileobj(file_in, file_out)



def get_freq_data(
        inTrunk,
        outTrunk,
        plinklocal
):
    """
    Wrapper for the --freq option in plink.
    """

    subprocess.run(
        [
            plinklocal,
            "--bfile", inTrunk,
            "--freq",
            "--out", outTrunk
        ],
        check = True
    )

    freq_file = str(outTrunk) + ".frq"

    # See if there are rare variants
    freq_data = pd.read_csv(
        freq_file,
        delimiter = '\s+'
    )

    return freq_data


def filter_excess_het(
        rule,
        markers,
        in_bedset,
        out_bedset,
        threshold = 0.1,
        sd = 2,
        result_file = '/dev/null',
        plot_file = False,
        plinklocal = None,
        rule_info = None
):
    """Runs plink het maf autosomal and produces output, including plot

    Wrapper around plink to do check sample heterozygosity and plot
    distribution.

    The function can be executed on two categories of markers, "common" and rare, that will trigger minimal and maximal maf thresholds, respectively.

    Saves results with saveYamlResults as well, and creates plots and various
    background files in the tmp area.

    """

    if markers == "common":
        maf = "--maf"
    elif markers == "rare":
        maf = "--max-maf"
    else:
        raise Exception(f"ERROR: Unexpected autosomal rarity {markers}.")
        return

    subprocess.run(
        [
            plinklocal,
            "--bfile", in_bedset,
            "--autosome",
            maf, str(threshold),
            "--het",
            "--out",
            out_bedset
        ],
        check = True
    )

    # We here have a .het file where low p-values for markers are to be removed
    het_p_values = out_bedset + ".het"
    compute_excess_het(het_p_values, out_bedset + ".exclude", sd)
    # compute_excess_het found .exclude for removal and .total for historgram

    subprocess.run(
        [
            plinklocal,
            "--bfile", in_bedset,
            "--remove", out_bedset + ".exclude",
            "--out", out_bedset,
            "--make-bed"
        ],
        check = True
    )

    dropouts = checkUpdates(
        in_bedset + ".fam",
        out_bedset + ".fam",
        cols = [0, 1],
        sanityCheck = "removal",
        fullList = True
    )

    dropouts.update(rule_info[rule])   # Metainfo and documentation about the rule
    dropouts["Threshold"] = f"{threshold} using distance {sd} standard deviations"
    dropouts["Rule"] = rule
    dropouts["Details"] = f"{out_bedset}.exclude with extension .details (filtered) and .total (all)"
    saveYamlResults(result_file, dropouts)
    title = (f'Sample heterozygosity ({markers} markes)\n'
             f'{dropouts.get("actionTakenCount")} outside threshold\n'
             f'--maf > {threshold} HET exceeds {sd} std.dev\n{dropouts.get("Timestamp")}')
    plot_hist(
        out_bedset + ".exclude.total",
        plot_file,
        column = "het_rate",
        title = title,
        separator = '\s+',
        threshold = 0,
        logx = False
    )
    return


def exclude_strand_ambigious_markers(
        input,
        output,
        plinklocal = None
):
    """ Runs plink to exlucde A/T and C/G SNPs
    input is the trunk of a .bed/.bim/.fam triplet
    ouput will be the corresponding set, without excluded markers
    The list over exluded markers can be found in ouput.excl
    plink will produce the output-bedset
    """
    try:
        df = pd.read_csv(
            input + ".bim",
            sep = '\s+',
            header = None
        ).astype(str)
    except Exception as e:
        print(f"Could not open .bim-file of bedset {input}, {str(e)}")
        return

    mask = ((df[4] == 'G') & (df[5] == 'C') |
            (df[4] == 'C') & (df[5] == 'G') |
            (df[4] == 'A') & (df[5] == 'T') |
            (df[4] == 'T') & (df[5] == 'A'))
    (df[mask])[1].to_csv(output + ".excl", index = False, header = False)

    subprocess.run(
        [
            plinklocal,
            "--bfile", input,
            "--exclude", output + ".excl",
            "--out", output,
            "--make-bed"
        ],
        check = True
    )



# def fix_rsid_map(mapfile, newmap):
#     """ NOT USED! Create a rsid mapping based on som Moba business-logic
#
#     Map of the rsid given in mapfile needs tweeking. newmap is produced
#     Multiple whitespaces in mapfile will become a single space in newmap
#     Could have been a lot more efficient
#
#     We do the following: (works for GSA/GSADM)
#     * Ignore lines that map to dot ('.')
#     * If multiple (comma-separated) ids are found in to-map, we use the first
#     * from strings containg .1 .2 ... .9 are ignored
#     """
#     try:
#         mappings = pd.read_csv(mapfile, usecols=[0, 1], names=['from', 'to'],
#                                delim_whitespace=True).astype(str)
#     except Exception as e:
#         print(f"Could not open file {mapfile}, {str(e)}")
#         return
#
#     with open(newmap, "w") as out:
#         multiAllele = re.compile("\.\d")
#         for index, row in mappings.iterrows():
#             # Elements to ignore
#             if row['to'] == ".":
#                 continue
#             if multiAllele.search(row['from']):
#                 continue   # eg rs222.1 is ignored
#             # Elements to simplify:
#             # 1 - Trucate everything including and after the first comma
#             to = re.sub(r",.+$", "", row['to'])
#             # Save
#             out.write(f"{row['from']} {to}\n")


def intersect_rsid(bim_small, bim_big, intersection, small_col = 1, big_col = 1):
    """Default assumes bim-ish files, that is tab-serarated columns.

    intersection is a file to be created
    Will create duplicate rsid if bim_big contains such
    If one of the files is large, pass that as bim_big for efficiency
    bim-files have rsid in 2. column (1, default). If your files
    contain rsid but in another format, pass the column number (0
    is first column)

    """
    my_name = inspect.currentframe().f_code.co_name      # this functions name
    m = 0   # max hits
    try:
        (smallDict, m) = dict_count_items(bim_small, [small_col], warn=False)
        with open(intersection, "w") as out:
            for line in open(bim_big):
                rsid = line.split()[big_col]
                if smallDict.get(rsid, 0) > 0:
                    out.write(f"{rsid}\n")
    except Exception as e:
        print(f"{my_name} Exception caught:  {str(e)}")


def copy_file(
        f,
        ext = ".back-up"
):
    """makes a .back-up file

    This is a debug-tool used keep a copy of a snakemake result that
    it would have deleted due to failure.

    """
    print(f"DEBUG> Making a backup of {f} to {f+ext}.")
    copyfile(f, f + ext)

def copy_bedset(
        trunkIn,
        trunkOut,
        extensions = [".bed", ".bim", ".fam"]
):
    """copies a plink bedset

    This takes a plink trunk as input and copies the files with the given extension

    """

    for extension in extensions:
        copyfile(trunkIn + extension, trunkOut + extension)

def dotplot(
        genomedata,
        prec = 2,
        x = 'x',
        y = 'y',
        c = 'c'
):
    """ Returns a plotnine object ready to be printet, but where extra lines can be added

    Assumes x and y are names of columns containing numbers, these will be rounded to precision decimals
    The preicision is there to not cluster set plot when there are two many points to be seen
    c is used for colouring and shaping - a colourblind palette will be used
    The data is found in a whitespace separated file where they x,y and c are headers

    """
    my_name = inspect.currentframe().f_code.co_name      # This function's name
    try:
        df = pd.read_csv(
            genomedata,
            usecols = [c, x, y],
            sep = '\s+'
        )
    except Exception as e:
        print(f"{my_name}: {str(e)}")
        return

    p = p9.ggplot(
        data = df.round(prec).drop_duplicates(),
        mapping = p9.aes(
            x = x,
            y = y,
            color = c,
            shape = c
        )
    )
    p += p9.scale_colour_brewer(
        type="qual",
        palette="Set1"
    )
    p += p9.geom_point()

    return p


def hweg_qq_plot(
        pfile,
        prec = 3,
        x = 'x'
):
    """Returns a plotnine object ready to be printed

    Values in column named x from file are picked up, and the plot
    will consist of -log(10) on the Y axis, and expected -log10 of a
    uniform distribution on the X axis

    The numbers will be rounded to prec(ision) and made unique before
    the plot is made.

    The precision is there to not cluster set plot when there are too
    many points to be seen The data is found in a whitespace separated
    pfile where they x denotes the header

    """
    my_name = inspect.currentframe().f_code.co_name      # This function's name
    try:
        df = pd.read_csv(
            pfile,
            usecols = [x],
            sep = '\s+'
        ).sort_values(
            by = [x],
            na_position = 'first',
            ascending = False
        )
    except Exception as e:
        print(f"{my_name}: {str(e)}")
        return
    # could/should have tested for 0 values here, avoiding log(0) problems
    df = -1 * np.log10(df)
    df.rename(
        columns = {x: "-log" + x},
        inplace = True
    )
    df['-logP_expected'] = -1 * np.log10(np.random.uniform(0, 1, len(df.index)))
    df['-logP_expected'] = df['-logP_expected'].sort_values(ascending = True).values

    p = p9.ggplot(
        data = df.round(prec).drop_duplicates(),
        mapping = p9.aes(
            y = '-logP',
            x = '-logP_expected'
        )
    )
    #p += p9.scale_colour_brewer(type="qual", palette="Set1")  # better for colourblind
    p += p9.geom_point()
    p += p9.geom_abline(
        slope = 1,
        intercept = 0,
        color = 'blue'
    )  # expexted hwe distibution

    return p


def sex_check(
        rule,
        in_bed,
        out_bed,
        f_threshold = 0.2,
        m_threshold = 0.8,
        config_sex_check_indep_pairwise = "20000 2000 0.5",
        result_file = '/dev/null',
        plot_file = False,
        plinklocal = None,
        rule_info = None
):
    """ Checks if sex according to .fam-file matches genotype. Removes mismatches

    Wrapper around several plink commands. in/out_bedset are plink .bed-files.

    Saves results with saveYamlResults as well, and creates plots and various
    backgroundfiles on the tmp-area.

    """

    tmpPath = Path(out_bed).parent
    inTrunk = plinkBase(in_bed)
    outTrunk = plinkBase(out_bed)

    # prune/chr23
    sex_check_parameters = config_sex_check_indep_pairwise.split()
    subprocess.run(
        [
            plinklocal,
            "--bfile", inTrunk,
            "--chr", "23",
            "--out", tmpPath/"pruned_sex_markers",
            "--indep-pairphase"
        ] + sex_check_parameters,
        check = True
    )

    subprocess.run(
        [
            plinklocal,
            "--bfile", inTrunk,
            "--extract", tmpPath/"pruned_sex_markers.prune.in",
            "--out", tmpPath/"pruned_for_sexcheck",
            "--make-bed"
        ]
    )

    # check sex on X chromosome
    subprocess.run(
        [
            plinklocal,
            "--bfile", tmpPath/"pruned_for_sexcheck",
            "--check-sex", str(f_threshold), str(m_threshold),
            "--out", tmpPath/"sexcheck_report_x",
            ],
        check = True
    )

    # find famid/id of the ones failing sexheck
    extract_list(
        tmpPath/"sexcheck_report_x.sexcheck",
        tmpPath/"sexcheck_report_x.remove",
        tmpPath/"sexcheck_report_x.details",
        colName = "^STATUS$",
        sep = None,
        condition = '==',
        threshold = "PROBLEM",
        key_cols = [0, 1],
        doc_cols = [0, 1, 5]
    )

    # remove these
    subprocess.run(
        [
            plinklocal,
            "--bfile", inTrunk,
            "--remove", tmpPath/"sexcheck_report_x.remove",
            "--out", outTrunk,
            "--make-bed"
        ],
        check=True
    )

    dropouts = checkUpdates(
        inTrunk + ".fam",
        outTrunk + ".fam",
        cols = [0, 1],
        sanityCheck = "removal",
        fullList = True
    )

    dropouts.update(rule_info[rule])   # Extra documentation about the rule
    dropouts["Threshold"] = f"Female={f_threshold} Male={m_threshold}"
    dropouts["Rule"] = rule
    saveYamlResults(result_file, dropouts)
    title = (f'Sex Check on chromosome X\n'
             f'{dropouts.get("actionTakenCount")} outside threshold\n'
             f'Threshold {dropouts.get("Threshold")}\n'
             f'{dropouts.get("Timestamp")}')
    plot_hist(
        tmpPath/"sexcheck_report_x.sexcheck",
        plot_file,
        column = "F",
        title = title,
        separator = '\s+',
        threshold = 0,
        logx = False,
        bins = 100
    )

    return


def egrep(pattern, in_file, out_file, switches=""):
    """ egrep wrapper. Lazy. Prone to path errors. It is pretty bad tbh.

        Consider simply using subprocess.run
    """
    subprocess.call(f'egrep {switches} {pattern} {in_file} > {out_file}',
                    shell=True)
    return


def count_families(famfile, regex):
    """Counts families with tuples in a fam-file (matching a certain pattern)

    Returns a dictionaly with key "nomatch" for the number of families
    in famfile with an id not matching regex.
    The keys 1,2,3 represents the number of matching families with 1,2
    or 3 members

    Writes a warning to stdout if families contain more than 3 memebers

    Use * as wildcard to check everything. This function is used to count
    the number of remaining families after we have tried to fix the pedigree

    """
    good_family = re.compile(regex)
    (all, n) = dict_count_items(famfile, cols=[0], warn=False)
    if n > 3:
        print("**** OUCH! We have families with up to {n} members!")
    counts = dict()   # counting occurences of 1,2,3 and hopefully not more
    for fam in all:
        if re.search(good_family, fam):  # For these we count tuples
            counts[all[fam]] = counts.get(all[fam], 0) + 1
        else:  # these do not match, but we count them still
            counts["nomatch"] = counts.get("nomatch", 0) + 1
    return(counts)

#
# def find_moba_pca_outlier(df):
#     """NOT IN USE/WORKING!
#
#     Actually this seems to work when memory is available.
#
#     The idea is to use stat_ellipse() with good parameters to draw an
#     ellipse/circle around the center - and remove them with at
#     corresponing pandas test df is a dataframe containing
#
#     * "PC1" and "PC2" columns (pca components, float values'
#     * "SuperPop" and "Population" columns. Will be edited for outliers
#       adding "(outlier")
#     A list of outliers will be produced to ... (file)
#
#     """
#     # print(df.head(40))
#     # need copy for this to work
#     df.loc[df["PC1"] > 0, 'SuperPop'] = df['SuperPop'] + "(outlier)"
#     df.loc[df["PC1"] > 0, 'Population'] = df['Population'] + "(outlier)"
#
#     p = p9.ggplot(
#         data = df,
#         mapping = p9.aes(
#             x = 'PC1',
#             y = 'PC2',
#             color = "Population"
#         )
#     )
#
#     # show an ellipse. sounds like a good idea, but didnt work for large sets
#     # while  harvest server was overloaded. If this works, we could use the
#     # corresponding test to remove outliers.
#     p += p9.stat_ellipse()
#     p += p9.geom_point()
#     p.save(
#         filename = "foo.png",
#         dpi = 300
#     )

# dir="/mnt/work2/gutorm/pipeOut/mod2-data-preparation/founders/"
# intersect_rsid("/mnt/work/gutorm/git/mobaGenetics-qc/qc-pipeline/snakefiles/foo", dir+"23", "bar", small_col=0, big_col=1)

def create_fam_map(fam_file, map_in_file,  map_out_file):
    """Creates a file that plink can use to rename individuals

    This is typically used to replace/obfuscate retrievalDetails_Ids
    with sentrixIds.  

    * fam_file is the original plink fam-file 

    * map_in_file is typically a subset of the samplesheet, no headers
      and only retrievalDetails_Ids and sentrixIds

    * map_out_file the mapping file expected later by plink --update-ids

    """
    # print(f"Doing: {fam_file} and {map_in_file} to {map_out_file}")
    my_name = inspect.currentframe().f_code.co_name  # Generic way of find function name
    try:
        # 4 columns 0-3
        fam = pd.read_csv(
            fam_file,
            header = None,
            sep = '\s+'
        )
        # 2 columns 4-5
        map = pd.read_csv(
            map_in_file,
            header = None,
            sep = '\s+'
        )
    except Exception as e:
        print(f"{my_name}: {str(e)}")
        return

    # print(map)
    all = fam.merge(
        map,
        left_on = [1],
        right_on = [0],
        indicator = True,
        validate = "1:1"
    )
    # 0_x is old familyname, 
    # 1_x is old ID (retrievalID), 1_y is new (sentrixId)
    all[["0_x","1_x","1_y","1_y"]].to_csv(
        map_out_file,
        sep = " ",
        header = False,
        index = False
    )
    # end create_fam_map



def restore_family_information(fam_files, batches, post_imputation_psam_file, new_psam_file):
    fam_dfs = []
    for i in range(len(fam_files)):
        fam_file = fam_files[i]
        df = pd.read_csv(fam_file, delim_whitespace=True, header=None, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])
        df["SID"] = batches[i]
        fam_dfs.append(df)
    fam_df = pd.concat(fam_dfs, ignore_index=True)
    psam_df = pd.read_csv(post_imputation_psam_file, delim_whitespace=True)
    psam_df = psam_df.reset_index()
    new_psam_df = pd.DataFrame(columns=["#FID", "IID", "SID", "PAT", "MAT", "SEX"])
    for index, psam_row in psam_df.iterrows():
        IID = psam_row["#IID"]
        IID_row = fam_df[fam_df["IID"] == IID]
        SID = IID_row["SID"].iloc[0]
        FID = IID_row["FID"].iloc[0]
        PAT = IID_row["PID"].iloc[0]
        MAT = IID_row["MID"].iloc[0]
        SEX = IID_row["Sex"].iloc[0]
        new_psam_df.iloc[index] = [FID, IID, SID, PAT, MAT, SEX]
    new_psam_df.to_csv(new_psam_file, sep="\t", index=False)


# merge is not implemented in plink 2 yet
# def merge_pgensets(pgens, out_trunk, plink2local):
#     """
#     pgen: list of .pgen-files
#     merges the pgen-sets associated with the .pgen-files in pgens into a single pgen-set with filebase out_trunk
#     """
#     cmd = f"""
#     # Generate list of files to merge
#     pgenset_dir=$(dirname "{out_trunk}")
#     echo {pgens} | tr ' ' '\\n' | sed 's/.pgen//' > $pgenset_dir/pgen_list.txt
#     {plink2local} --pmerge-list $pgenset_dir/pgen_list.txt --out {out_trunk}
#     """
#     subprocess.run(cmd, shell=True, check=True)


def main():
    # if you want to test a function
    print("Main called")
    y = checkUpdates(
        "/mnt/work/gutorm/1000Genomes/all_phase3.bim",
        "/mnt/work2/gutorm/pipeOut/mod3-good-markers/pca_ref.bim",
        cols = [0, 1, 3, 4, 5],
        indx = 1,
        sanityCheck = "none",
        fullList = True,
        allele_flip = True
    )
    print(y)
if __name__ == "__main__":
    main()
