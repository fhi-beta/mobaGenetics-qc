#!/usr/bin/env python3
import pandas as pd
import yaml
from pathlib import Path
import sys
import argparse
import re

def create_new_flag_col(df, c, qcFile,  failPass=["Failed","Passed"]):
    """ Append a row to the flag matrix

    df is the frame to be augmentet by a new columne named c
    the results from the QC are samples found in qcFile
    When a match is not found, the corresponding string from failPass is used
    qcFile is typical result from checkUpdates() with fullList=True. If the sample/marker
    has been modified, it is present in the "xitems" list of the .details file. 

    Dropped samples should be marked by a X. When that happens, all following columns will
    be set to x (lower case). 
    
    """
    try:    # creating a dictionary of samples/markers listed
        with open(qcFile) as file:
            qcFull = yaml.load(file, Loader=yaml.FullLoader)
        allItems = qcFull["xitems"]
        item = re.compile('^\w+')  # sample or marker id
        items = [item.match(i).group() for i in allItems]
        qc = dict.fromkeys(items,1)
    except Exception as e:
        print(f"Could not create dictionary from yaml file {qcFile} ({str(e)})")
        return

    # we translate the qc list of samples to a pass/not passed code
    results = list(df["Id"])             # These are the samples/markers
    for i, result in enumerate(results): # found in the qc-fresults or not?
        results[i] = failPass[0] if qc.get(result,0) == 0 else failPass[1]
    df[c] = results                      # add a new column
    # If the sample was previously dropped, the second last column will cointain x or X.
    # In that case, we just add an 'x': Not being in the .detail file means nothing.
    if len(df.columns) > 1:
        df.loc[df.iloc[:,-2].str.upper() == 'X', [c]] = 'x'


def main(argv):
    """
    Creates a flaglist, based on an exiting samplefile.
    The filesamplefile should sport names used through the whole qc process - typicall after having lab-ids updated to moba ids.
    Is dependent on a configfile (-c) that shows what columns to make and call them, as well as results from the qc-pipeline.
    A default configfile flagTemplate.csv is provided at the installdirecotry. 
    

    The program could be used for markers as well, but that has not been tried
    Displays a usage if run without parameters.

    """
    
    parser = argparse.ArgumentParser(description='Creates a flag-list showing all QC results for all samples')
    parser.add_argument("--samplefile","-s", required=True, help="Typically a plink .fam file. parents-updatet.fam is usually a good one as it is early in the pipe")
    parser.add_argument("--configfile","-c", default="flagTemplate.csv", help="A csv-file describing columns to add to the flagfile")
    parser.add_argument("--resultdir","-r", required=True, help="A direcory containing QC results")
    parser.add_argument("--output","-o", required=False, help="Result for the flag-matrix. Default sampleFlags.txt on resuldir")
    parser.add_argument("--verbatim","-v", default=0, type=int, required=False, help="Chatty during execution. Default 0, any other number is true")
    args = parser.parse_args()
    chatty = args.verbatim != 0
    result_file = args.output
    if result_file == None:
        result_file = Path(args.resultdir)/"sampleFlags.txt"

    try:
        # New columns to add, as well as their values
        names = pd.read_csv(args.configfile, sep=":", comment="#", skipinitialspace=True)
        # The original sample list
        all = pd.read_csv(args.samplefile, header=None, sep=" ", 
                          usecols=[0,1], names=["Fam","Id"])
    except Exception as e:
        print(f"Sorry could not do it. {str(e)}")

    last_mod_time = 0
    last_file = "No such file"
    
    for col_name in names.columns:    # each iteration adds a new column to the flag_matrix
        if chatty : print (col_name)
        try:
            qc = args.resultdir / Path(names.loc[0,col_name])
        except Exception as e:
            print(f"Error: Failing to read resultfile for {col_name}: {names.loc[0,col_name]}. {str(e)}")
            print (names)
        if qc.stat().st_mtime < last_mod_time:
            sys.stderr.write(f"Warning {qc.name} is older than {last_file}\n")
        create_new_flag_col(all, col_name, qcFile=qc, 
                   failPass=[names.loc[2,col_name],names.loc[1,col_name]])  # These are pass/fail strings
        last_mod_time = qc.stat().st_mtime
        last_file = qc.name

    all.to_csv(result_file, sep = ' ', index=False)

if __name__ == "__main__":
   main(sys.argv[1:])
