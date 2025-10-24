# preparing for phasing check (work in progress)

import pysam
from datetime import datetime


chr = "22"
best_snps = "/mnt/archive3/phasing_test/best_snps"
expected_relationships_file = "/mnt/archive3/phasing_test/shapeit_fam"
bcf_file = "/mnt/archive3/phasing_test/chr22.ac.bcf"
children = []
family_relations = {}
n = 500

variant_positions = []
with open(best_snps, 'r') as f:
    for line in f:
        chrom, pos, id = line.strip().split('\t')
        pos = int(pos)
        if chrom == chr:
            variant_positions.append((chrom, pos, id))

with open(expected_relationships_file) as f:
    for line in f:
        child, father, mother = line.strip().split()
        if child != "NA" and father != "NA" and mother != "NA":
            children.append(child)
            family_relations[child] = {'father': father, 'mother': mother}

children_data = {}
for child in children:
    children_data[child] = {"e":0, "n":0}

vcf_in = pysam.VariantFile(bcf_file)

def prune_list(original_list, n):
    """
    Prunes the original list to a size of n by removing elements equally spaced apart.
    """
    if n > len(original_list):
        return original_list
    step = len(original_list) / n
    pruned_list = [original_list[int(i * step)] for i in range(n)]
    return pruned_list

variant_positions_pruned = prune_list(variant_positions, n)



for variant in variant_positions_pruned:
    chrom = variant[0]
    pos = variant[1]
    id = variant[2]
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    for record in vcf_in.fetch(chrom, pos - 1, pos):
        for child_id in children:
            child_genotype = record.samples.get(child_id, {}).get('GT')
            if (child_genotype in [(0, 1), (1, 0)]):
                father_id = family_relations[child_id]['father']
                father_genotype = record.samples.get(father_id, {}).get('GT')
                if father_genotype in [(0, 0), (1, 1)]:
                    mother_id = family_relations[child_id]['mother']          
                    mother_genotype = record.samples.get(mother_id, {}).get('GT')
                    if mother_genotype in [(0, 0), (1, 1)]:
                        if father_genotype != mother_genotype:
                            children_data[child_id]['n'] += 1
                            if (child_genotype[0] == 0 and father_genotype != (0,0)) or (child_genotype[0] == 1 and father_genotype != (1,1)):
                                children_data[child_id]['e'] += 1

output_file_path = "/mnt/archive3/phasing_test/phasing_check_test2" 
with open(output_file_path, 'w') as output_file:
    for child in children:
        n = children_data[child]['n']
        e = children_data[child]['e']
        if n > 0:
            r = e/n
        else:
            r = "NA"
        output_line = f"{child}\t{e}\n{n}\tr"
                    