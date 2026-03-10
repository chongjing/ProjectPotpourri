import sys
import csv

def reformat_go(input_file, output_file):
    protein_go_map = {}

    # Read the input file
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            protein = row[0]
            go_terms = row[1].split(',')
            if protein not in protein_go_map:
                protein_go_map[protein] = set()
            protein_go_map[protein].update(go_terms)

    # Write the output file
    with open(output_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for protein, go_terms in protein_go_map.items():
            writer.writerow([protein, ','.join(sorted(go_terms))])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 reformat_GO.py input.tsv output.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reformat_go(input_file, output_file)
