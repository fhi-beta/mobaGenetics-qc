# script for checking phasing and Mendelian errors for trios

import argparse
import pysam
from datetime import datetime



def check_all(record, f):
    f['n_missing'] += 1
    child_id = f['child']
    child_genotype = record.samples.get(child_id, {}).get('GT')
    father_id = f['father']
    father_genotype = record.samples.get(father_id, {}).get('GT')
    mother_id = f['mother']          
    mother_genotype = record.samples.get(mother_id, {}).get('GT')
    proceed = True

    if child_genotype is None or child_genotype[0] is None or child_genotype[1] is None:
        f['e_child_missing'] += 1
        proceed = False
    if father_genotype is None or father_genotype[0] is None or father_genotype[1] is None:
        f['e_father_missing'] += 1
        proceed = False
    if mother_genotype is None or mother_genotype[0] is None or mother_genotype[1] is None:
        f['e_mother_missing'] += 1
        proceed = False

    if proceed:
        f['n_mendel'] += 1  
        if not ((child_genotype[0] in father_genotype and child_genotype[1] in mother_genotype) or (child_genotype[1] in father_genotype and child_genotype[0] in mother_genotype)):
            f['e_mendel'] += 1
        else:
            if child_genotype in [(0, 1), (1, 0)]:
                if (father_genotype in [(0, 0), (1, 1)] or mother_genotype in [(0, 0), (1, 1)]):
                    f['n_phasing'] += 1
                    if (child_genotype[0] == 0 and father_genotype == (1, 1)) or (child_genotype[0] == 1 and father_genotype == (0, 0)) or \
                       (child_genotype[1] == 0 and mother_genotype == (1, 1)) or (child_genotype[1] == 1 and mother_genotype == (0, 0)):
                        f['e_phasing'] += 1
                if (father_genotype in [(0, 0), (1, 1)] and mother_genotype in [(0, 0), (1, 1)] and mother_genotype != father_genotype):
                    f['n_phasing_hom'] += 1
                    if (child_genotype[0] == 0 and father_genotype != (0, 0)) or (child_genotype[0] == 1 and father_genotype != (1, 1)):
                        f['e_phasing_hom'] += 1
    
def check_phasing_hom(record, f):
    child_id = f['child']
    child_genotype = record.samples.get(child_id, {}).get('GT')
    if child_genotype is not None and child_genotype[0] is not None and child_genotype[1] is not None:
        if child_genotype in [(0, 1), (1, 0)]:
            father_id = f['father']
            father_genotype = record.samples.get(father_id, {}).get('GT')
            if father_genotype is not None and father_genotype[0] is not None and father_genotype[1] is not None:
                mother_id = f['mother']          
                mother_genotype = record.samples.get(mother_id, {}).get('GT')
                if mother_genotype is not None and mother_genotype[0] is not None and mother_genotype[1] is not None:
                    if (father_genotype in [(0, 0), (1, 1)] and mother_genotype in [(0, 0), (1, 1)] and mother_genotype != father_genotype):
                        f['n_phasing_hom'] += 1
                        if (child_genotype[0] == 0 and father_genotype != (0, 0)) or (child_genotype[0] == 1 and father_genotype != (1, 1)):
                            f['e_phasing_hom'] += 1

def main(args):
    trios_file = args.trios
    bcf_file = args.bcf
    output = args.output
    optimize = args.optimize
    pf = args.pf
    header = args.header
    trios = []
    with open(trios_file) as f:
        if header:
            next(f)
        for line in f:
            child, father, mother, iid_batch, pat_batch, mat_batch, iid_chip, pat_chip, mat_chip, iid_reg, pat_reg, mat_reg, parents_in_batch, shared_chips, parents, orig_batch, orig_parents_in_batch, move_from, move_to, moved= line.strip().split()
            if child != "NA" and father != "NA" and mother != "NA":
                trios.append({'child':child, 'father': father, 'mother': mother, "reg_id": iid_reg, "parents_in_batch":parents_in_batch,"orig_parents_in_batch":orig_parents_in_batch, "shared_chips": shared_chips, "e_phasing":0, "n_phasing":0, "e_phasing_hom":0, "n_phasing_hom":0, "e_child_missing":0, "e_father_missing":0, "e_mother_missing":0, "n_missing":0, "e_mendel":0, "n_mendel":0})

    vcf_in = pysam.VariantFile(bcf_file)
    counter = 0
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(current_time)
    for record in vcf_in:
        if counter % pf == 0:
            id = record.id
            chrom = record.chrom
            pos = record.pos
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print(f"Time: {current_time} Count: {counter} Chrom: {chrom}, Pos: {pos}, ID: {id}")
        counter+=1
        if optimize:
            for f in trios:
                check_phasing_hom(record, f)
        else:
            for f in trios:
                check_all(record, f)

    with open(output, 'w') as output_file:
        if optimize:
            output_line = "iid\tpat\tmat\treg_id\tparents_in_batch\torig_parents_in_batch\tshared_chips\te_phasing_hom\tn_phasing_hom\tr_phasing_hom\n"
        else:
            output_line = "iid\tpat\tmat\treg_id\tparents_in_batch\torig_parents_in_batch\tshared_chips\te_phasing\tn_phasing\tr_phasing\te_phasing_hom\tn_phasing_hom\tr_phasing_hom\te_mendel\tn_mendel\tr_mendel\te_child_missing\te_father_missing\te_mother_missing\tn_missing\tr_child_missing\tr_father_missing\tr_mother_missing\n"
        output_file.write(output_line)
        for f in trios:
            child = f['child']
            father = f['father']
            mother = f['mother']
            reg_id = f['reg_id']
            shared_chips = f['shared_chips']
            parents_in_batch = f['parents_in_batch']
            orig_parents_in_batch = f['orig_parents_in_batch']
            e_phasing_hom = f['e_phasing_hom']
            n_phasing_hom = f['n_phasing_hom']
            if n_phasing_hom > 0:
                r_phasing_hom = e_phasing_hom/n_phasing_hom
            else:
                r_phasing_hom = "NA"
            if not optimize:
                e_phasing = f['e_phasing']
                n_phasing = f['n_phasing']
                if n_phasing > 0:
                    r_phasing = e_phasing/n_phasing
                else:
                    r_phasing = "NA"
                e_child_missing = f['e_child_missing']
                e_father_missing = f['e_father_missing']
                e_mother_missing = f['e_mother_missing']
                n_missing = f['n_missing']
                e_mendel = f['e_mendel']
                n_mendel = f['n_mendel']
                if n_mendel > 0:
                    r_mendel = e_mendel/n_mendel
                else:
                    r_mendel = "NA"
                if n_missing > 0:
                    r_child_missing = e_child_missing/n_missing
                    r_father_missing = e_father_missing/n_missing
                    r_mother_missing = e_mother_missing/n_missing
                else:
                    r_child_missing = "NA"
                    r_father_missing = "NA"
                    r_mother_missing = "NA"
            if optimize:
                output_line = f"{child}\t{father}\t{mother}\t{reg_id}\t{parents_in_batch}\t{orig_parents_in_batch}\t{shared_chips}\t{e_phasing_hom}\t{n_phasing_hom}\t{r_phasing_hom}\n"
            else:
                output_line = f"{child}\t{father}\t{mother}\t{reg_id}\t{parents_in_batch}\t{orig_parents_in_batch}\t{shared_chips}\t{e_phasing}\t{n_phasing}\t{r_phasing}\t{e_phasing_hom}\t{n_phasing_hom}\t{r_phasing_hom}\t{e_mendel}\t{n_mendel}\t{r_mendel}\t{e_child_missing}\t{e_father_missing}\t{e_mother_missing}\t{n_missing}\t{r_child_missing}\t{r_father_missing}\t{r_mother_missing}\n"
            output_file.write(output_line)

    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(current_time)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Phasing check')
    parser.add_argument('--trios', type=str, required=True, help='Path to trios file')
    parser.add_argument('--bcf', type=str, required=True, help='Path to BCF file')
    parser.add_argument('--output', type=str, required=True, help='Output file path')
    parser.add_argument('--optimize', action='store_true', help='Check only phasing')
    parser.add_argument('--header', action='store_true', help='Skip header')
    parser.add_argument('--pf', type=int, required=False, default=1000, help='Print frequency (defaults to 1000)')
    args = parser.parse_args()
    main(args)