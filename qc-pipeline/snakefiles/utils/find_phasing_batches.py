import argparse
def main(args):
    debug = False
    if debug:
        families_file = "/mnt/archive3/phasing_test/debug/families"
        n_batches = 2
        output_trunk = "/mnt/archive3/phasing_test/debug/phasing_batch"
    else:
        families_file = args.families
        n_batches = args.n_batches
        output_trunk = args.output_trunk

    families = []
    with open(families_file, 'r') as file:
        for line in file:
            family_members = line.strip().split("\t")
            families.append(family_members)

    phasing_batches = [[] for i in range(n_batches)]
    for i, family in enumerate(families):
        batch_index = i % n_batches
        phasing_batches[batch_index].extend(family)

    for i in range(n_batches):
        filepath = f"{output_trunk}.batch{i}"
        samples = phasing_batches[i]         
        with open(filepath, 'w') as file:
            for sample in samples:
                file.write(f"{sample}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find phasing batches')
    parser.add_argument('--families', type=str, required=True, help='Path to input families file')
    parser.add_argument('--n_batches', type=int, required=True, help='Number of phasing batches to create')
    parser.add_argument('--output_trunk', type=str, required=True, help='Path prefix for output phasing batch files')
    args = parser.parse_args()
    main(args)