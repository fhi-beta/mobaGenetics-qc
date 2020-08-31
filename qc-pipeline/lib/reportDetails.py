#!/usr/bin/env python3
import pandas as pd
import yaml
from pathlib import Path
from subprocess import PIPE, run
import sys
import argparse
import re


# def create_new_flag_col(df, c, qcFile, failPass=["Failed", "Passed"]):
#     """ Append a row to the flag matrix

#     df is the frame to be augmentet by a new columne named c
#     the results from the QC are samples found in qcFile
#     When a match is not found, the corresponding string from failPass is used
#     qcFile is typical result from checkUpdates() with fullList=True. If the sample/marker
#     has been modified, it is present in the "xitems" list of the .details file.

#     Dropped samples should be marked by a X. When that happens, 
#     all following columns will
#     be set to x (lower case).

#     """
#     try:    # creating a dictionary of samples/markers listed
#         with open(qcFile) as file:
#             qcFull = yaml.load(file, Loader=yaml.FullLoader)
#         allItems = qcFull["xitems"]
#         item = re.compile('^\w+')  # sample or marker id
#         items = [item.match(i).group() for i in allItems]
#         qc = dict.fromkeys(items, 1)
#     except Exception as e:
#         print(f"Could not create dictionary from yaml file {qcFile} ({str(e)})")
#         return
#     # we translate the qc list of samples to a pass/not passed code
#     results = list(df["Id"])              # These are the samples/markers
#     for i, result in enumerate(results):  # found in the qc-fresults or not?
#         results[i] = failPass[0] if qc.get(result, 0) == 0 else failPass[1]
#     df[c] = results                      # add a new column
#     # If the sample was previously dropped, the second last column contains x or X.
#     # In that case, we just add an 'x': Not being in the .detail file means nothing.
#     if len(df.columns) > 1:
#         df.loc[df.iloc[:, -2].str.upper() == 'X', [c]] = 'x'

def parse_yaml(file, value):
    """ Checking a result-file for a marker/sample

    file is the file to check, and is assumed to have been made by the
    qc-pipleine and contains a yaml-structure.

    If the parameter value is found (using egrep), the file will be parsed 
    and selcted items from the yaml-file will be printed out.
    This is values as the rule that has processed the marker, when it was run
    and so on.
    
    The yaml-file is expected to contain the marker/sample only if the 
    rule removed/replaced it.

    Returns the number of matches found, -1 on eggrep error
    """
    #print (value, file)
    # Counting matches before we start parsing the file
    result = run(["egrep", "-c", value, file], stdout=PIPE, stderr=PIPE,
                 universal_newlines=True)
    if result.returncode > 1:
        print("egrep returned error")
        return -1
    elif result.returncode == 1: # no match
        return 0
    # We here have a match. Parse the yaml-file and dump results.
    with open(file, 'r') as stream:
        r = yaml.safe_load(stream)  # r for result

    print (f'Item {r["rule action"]}d by rule {r["Rule"]} ({r["Rule order"]})'
           f' on {r["Timestamp"]}')
    print (f'Test/description: {r["QC test"]}\n{r["Description"]}')
    print (f'See {file} for full details')
    

    # end parse_yaml
def main(argv):
    """Reports on a sample or marker

    Displays a usage if run without parameters.

    """

    result_dir = "."
    files_to_check = "markers.yaml"
    parser = argparse.ArgumentParser(description='Reports on what rules in the QC did what to a sample/marker')
    parser.add_argument("--resultdir","-r", default=result_dir, help=f"A direcory containing QC results (default {result_dir})")
    parser.add_argument("--configfile","-c", default=files_to_check, help=f"A list of rules result files to check (default {files_to_check}). Order of files should match expected pipleine age, if not use the --silent switch")
    parser.add_argument("--id","-i", required=True, help=f"identificator to look for in the files supplied by -c. For markers, you can choose between snpid or 'chr:pos:a1:a2' ")

    parser.add_argument("--verbatim","-v", default=0, type=int, required=False,help="Chatty during execution. Default 0, any other number is true")
    parser.add_argument("--silent","-s", default=0, type=int, required=False,help="Silent. Default 0, displays warning. Set to 1 to supress warnings")
    args = parser.parse_args()
    chatty = args.verbatim != 0
    silent = args.silent != 0
    
    print("Should be able to check that the sample/markers is in the final set...")
    target = args.id
    result_dir = args.resultdir
    files_to_check = args.configfile
    
    is_marker = target.split(":")
    if  len(is_marker)  == 4 :
        # This is a chr:pos:a1:a1 syntax later egrep will have probles
        if not silent:
            print ("Looking for a postion based marker")
        # Building a egrep regexp matching .bim files. Quote the \S (non-blank)
        target = target.split(":")[0] + " \\S+ " + target.split(":")[1] + \
                 " " + target.split(":")[2] + " " + target.split(":")[3]
            
    if chatty:
        print(f"Looking for '{target}' in results found in directory '{result_dir}' "
              f"Checking only files listed in {files_to_check}"  )

    try:
        # Grab files to check|
        res_files = pd.read_csv(args.configfile, sep=":", comment="#", skipinitialspace=True,names=["result_file"] )
    except Exception as e:
        print(f"Sorry could not do it. {str(e)}")


        
    last_mod_time = 0
    last_file = "No such file"
    # print (res_files)
    for index, row in res_files.iterrows():  # process all files
        res_file_name = row["result_file"]
        if chatty:
            print(f"Checking {res_file_name}")
        try:
            qc_output = args.resultdir / Path(res_file_name)
        except Exception as e:
            print(f"Error: Failing to read resultfile for {qc_output}: {str(e)}")
        if qc_output.stat().st_mtime < last_mod_time and not silent:
            sys.stderr.write(f"Warning {res_file_name} is older than {last_file}\n")
            
        parse_yaml(qc_output, target)
        last_mod_time = qc_output.stat().st_mtime
        last_file = qc_output

    
if __name__ == "__main__":
    main(sys.argv[1:])
