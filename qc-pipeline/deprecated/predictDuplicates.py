# Primitive duplicate check. See 'docs' under predict_and_list_duplicates

import pandas as pd
import yaml

# samples mapping
recode =  "/mnt/archive/ROTTERDAM2/delivery-fhi/data/auxiliaries/recode-files/recode-files-original-fhi/recode-ids-rotterdam2.txt"
#marker mapping
recode = "/mnt/archive/ROTTERDAM2/delivery-fhi/data/auxiliaries/recode-files/recode-rsid.txt"

#  datafile
data = "/mnt/work/gutorm/qcTest/qcrot2/fullNewOutput/mod1-data-preparation/tmp/dup_markers_removed.bim"
data = "/mnt/work/gutorm/qcTest/qcrot2/fullNewOutput/mod1-data-preparation/tmp/aa_theta_removal.bim"

resultfile = "dups.yaml"

def dict_head(d, n):
    """
    Returns a subset of dict, n elements. Mostly for dubugging/pretty printing.
    """
    return(dict(list(d.items())[:n-1]))


def predict_and_list_duplicates(mapping_file, data_file,
                                mapping_filecols=[0,1], data_filecol=1):
    """
    Assume that mapping_file is used on data_file: Would that create duplicate ids and wich ones?
    mapping_cols are the from and to columns from mapping_file
    data_file(column data_filecol) is only checked if mapping to/from would create duplicates:
    That is if items separate elements (sample/marker-id) would be merged.
    The mappinglist in itself could theoretically create duplicates without this being
    a problem for data_file.
    returns a dictionary of the problematic samples/markers: Key is what the values in data_file will be mapped to
    """

    try: 
        df = pd.read_csv(mapping_file, delim_whitespace=True, header=None,
                         names=["from","to"], usecols=mapping_filecols).astype(str)
    except Exception as e:
        print(f"Could not open mapping file {mapping_file}, {str(e)}")

    to_from =  df.set_index('to').to_dict()['from']
    from_to =  df.set_index('from').to_dict()['to']
    if len(from_to) == len(to_from):
        print(f"Unique items in both in/out set:  {len(from_to)}")
        return
    print(f"Mapping in {mapping_file} is not one-to-one: Unique items from/to: {len(from_to)}/{len(to_from)}")

    # Now we have to find sample/marker-id that will create crash/lose sample-info

    try:
        data = pd.read_csv(data_file, usecols=[data_filecol], squeeze=True,
                       delim_whitespace=True, header=None)
    except Exception as e:
        print(f"Could not open data file {data_file}, {str(e)}")       

    dupDict = dict()  # Dictionary of duplicates (and non-duplicates)
    nonMatches = 0
    for i in data:
        map_to = from_to.get(i,"Argh. No match")
        if map_to == "Argh. No match":
            nonMatches +=1
            continue
        if map_to in dupDict : 
            if i not in dupDict[map_to]:  # Remember i for later, but only once
                dupDict[map_to].append(i)
        else: # Not seen before
            dupDict[map_to] = [i]

    if nonMatches > 0 : print(f"{nonMatches} items from {data_file} will not be mapped/translated when using {mapping_file}")
    for i in list(dupDict.keys()):   # Clean up non-duplicates
        if len(dupDict[i]) == 1 : del dupDict[i]

    return dupDict
        
print(f"Prediciting mapping of {data} by using the mapping {recode}")    
problems = predict_and_list_duplicates(recode, data)
print(f"The mapping will create {len(problems)} duplicates, shown in {resultfile}. \n"
      f"  foo:\n"
      f"  - oldfoo\n"
      f"  - oldfoo2\n"
      f"in {resultfile} means that ids oldfoo and oldfoo2 from {data} both will become foo after applying {recode}"
)

yaml.dump(problems, open(resultfile,"w"))
