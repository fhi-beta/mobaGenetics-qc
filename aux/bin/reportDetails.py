#!/usr/bin/env python3
import pandas as pd
import yaml
from pathlib import Path
from subprocess import PIPE, run
import sys
import argparse
import re



def parse_yaml(file, value):
    """Checking a result-file for a marker/sample

    file is the file to check, and is assumed to have been made by the
    qc-pipleine and contains a yaml-structure.

    If the parameter value (treated as regexp) is found (using egrep), the file will be
    parsed and selcted items from the yaml-file will be printed out.
    These are details like values as the name of the rule that has
    processed the marker, when it was run and so on.
    
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
    match_details = ""
    if int(result.stdout) > 1:
        match_details = f"!{int(result.stdout)} matches found!"
    with open(file, 'r') as stream:
        r = yaml.safe_load(stream)  # r for result

    print (f'Item {r["rule action"]}d by rule {r["Rule"]} ({r["Rule order"]})'
           f' on {r["Timestamp"]} {match_details}')
    print (f'Test/description: {r["QC test"]}\n{r["Description"]}')
    check_for = "Callrates"
    if check_for in r:
        print (f'{check_for}: {r[check_for]}')
    print (f'See {file} for full details\n')
    

    # end parse_yaml
def main(argv):
    """Reports on a sample or marker

    Displays a usage if run without parameters.

    """

    result_dir = "."
    files_to_check = "markers.yaml"
    parser = argparse.ArgumentParser(description='Reports on what rules in the QC did what to a sample/marker')
    parser.add_argument("--resultdir","-r", default=result_dir, help=f"A direcory containing QC results (default {result_dir}). Prefix to all files listet in configfile")
    parser.add_argument("--configfile","-c", default=files_to_check, help=f"A list of rules result files to check (default {files_to_check}). Order of files should match expected pipleine age, if not use the --silent switch")
    parser.add_argument("--id","-i", required=True, help=f"identificator to look for in the files supplied by -c. For markers, you can choose between snpid or 'chr:pos:a1:a2' ")

    parser.add_argument("--verbatim","-v", default=0, type=int, required=False,help="Chatty during execution. Default 0, any other number is true")
    parser.add_argument("--silent","-s", default=0, type=int, required=False,help="Silent. Default 0, displays warning. Set to 1 to supress warnings")
    args = parser.parse_args()
    chatty = args.verbatim != 0
    silent = args.silent != 0
    
    print("Should be able to check that the sample/markers is in the final set...")
    target = args.id
    result_dir = Path(args.resultdir)
    files_to_check = args.configfile
    
    is_marker = target.split(":")
    if  len(is_marker)  == 4 :
        # This is a chr:pos:a1:a1 syntax later egrep will have problems
        if not silent:
            print ("Looking for a postion based marker")
        # Building an egrep regexp matching .bim files. Quote the \S (non-blank)
        target = target.split(":")[0] + " \\S+ " + target.split(":")[1] + "("\
            # The following two lines (note the or |) handle allele flips
                 "( " + target.split(":")[2] + " " + target.split(":")[3] + ")|" +\
                 "( " + target.split(":")[3] + " " + target.split(":")[2] + ")" +\
                 ")"
    if chatty:
        print(f"Looking for '{target}'  "
              f"Checking only files listed in {files_to_check} (prefixed with '{result_dir}')"  )

    try:
        # Grab files to check|
        res_files = pd.read_csv(args.configfile, sep=":", comment="#",
                                skipinitialspace=True,names=["result_file"] )
    except Exception as e:
        print(f"Sorry could not do it. {str(e)}")


        
    last_mod_time = 0
    last_file = "No such file"
    # print (res_files)
    for index, row in res_files.iterrows():  # process all files
        res_file_name = result_dir / Path(row["result_file"])
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
