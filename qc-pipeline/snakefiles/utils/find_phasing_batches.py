# Script for splitting the samples into n_phasing_batches approximaely equal batches where all related samples are in the same batch

import pandas as pd
import numpy as np
from datetime import datetime
import copy
rel_file = "/mnt/archive3/phasing_test/phase_merged_reshuffle/updated_relations"
output_trunk =  "/mnt/archive3/phasing_test/phasing_batch"
fid_file = "/mnt/archive3/phasing_test/fids"
n_phasing_batches = 3
updated_rel = pd.read_csv(rel_file, sep = "\t", low_memory=False)


all_children = updated_rel[updated_rel['parents'] > 0]['iid'].tolist()
all_samples = updated_rel['iid'].tolist()
delegated_children = []
all_fathers = updated_rel['pat'].dropna().unique().tolist()
delegated_fathers = []
all_mothers = updated_rel['mat'].dropna().unique().tolist()
all_parents = all_fathers + all_mothers
delegated_mothers = []
not_children = list(set(all_samples) - set(all_children))
unrelated = list(set(not_children) - set(all_parents))
current_batch = 0
phasing_batches = [[] for i in range(n_phasing_batches)]

counter = 0
pf = 1000

def split_list(lst, n):
    k, m = divmod(len(lst), n)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]


def extend_no_duplicates(list1, list2):
    set1 = set(list1)
    list1.extend(item for item in list2 if item not in set1)
    return list1


def find_parents_and_siblings(child, previously_found = [[],[],[]]):
    total_found = copy.deepcopy(previously_found)
    mother = updated_rel.loc[updated_rel['iid'] == child, 'mat'].values[0]
    father = updated_rel.loc[updated_rel['iid'] == child, 'pat'].values[0]
    if not (isinstance(mother, float) and np.isnan(mother)):
        mother_offspring = updated_rel.loc[updated_rel['mat'] == mother, 'iid'].tolist()
        mother_list = [mother]
        delegated_mothers.append(mother)
    else:
        mother_offspring = []
        mother_list =[]
    if not (isinstance(father, float) and np.isnan(father)):
        father_offspring = updated_rel.loc[updated_rel['pat'] == father, 'iid'].tolist()
        father_list = [father]
        delegated_fathers.append(father)
    else:
        father_offspring = []
        father_list = []
    siblings = father_offspring
    siblings.extend(mother_offspring)
    siblings = list(set(siblings))
    total_found[0] = extend_no_duplicates(total_found[0], siblings)
    total_found[1] = extend_no_duplicates(total_found[1], father_list)
    total_found[2] = extend_no_duplicates(total_found[2], mother_list)
    for sibling in siblings:
        if sibling != child and (not sibling in previously_found[0]):
            total_found = find_parents_and_siblings(sibling, total_found)
    return total_found

phasing_batches = split_list(unrelated, n_phasing_batches)
families = [[un] for un in unrelated]
for sample in all_children:
    if counter % pf == 0:
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        n_delegated_samples = sum([len(batch) for batch in phasing_batches])
        print(f"{current_time}: Processing child {counter}. Total delegated samples: {n_delegated_samples}")
    counter += 1
    if sample not in delegated_children:
        child = sample
        family_lists = find_parents_and_siblings(child)
        children = family_lists[0]
        father_list = family_lists[1]
        mother_list = family_lists[2]
        family = children
        family.extend(father_list)
        family.extend(mother_list)
        families.append(family)
        delegated_children.extend(children)
        phasing_batches[current_batch] = phasing_batches[current_batch] + family
        if current_batch == n_phasing_batches -1:
            current_batch = 0
        else:
            current_batch += 1

fid_trunk = "fid_"

with open(fid_file, 'w') as file:
    for i in range(len(families)):
        family = families[i]
        fid = f"{fid_trunk}{i}"
        for sample in family:
            line = f"{fid}\t{sample}\n"
            file.write(line) 

for i in range(n_phasing_batches):
    filepath = f"{output_trunk}.batch{i}"
    samples = phasing_batches[i]         
    with open(filepath, 'w') as file:
        for sample in samples:
            file.write(f"{sample}\n")