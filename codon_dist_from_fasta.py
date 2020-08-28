import os
import argparse
import codon_dist as cd


def parse_fasta(file_handle):
    uids = []
    cds_list = []
    for line in file_handle:
        if line.startswith('>'):
            uid = line.split(';')[0][5:]
            uids.append(uid)

        else:
            seq = line.strip()
            cds_list.append(seq)

    return uids, cds_list


parser = argparse.ArgumentParser()
parser.add_argument("infile", type=argparse.FileType('r'), help="fasta file from get_species_info_by_uid.py")
parser.add_argument("target_directory", help="directory to store output files")
output = parser.add_mutually_exclusive_group(required=True)
output.add_argument("-o", "--outfile", help="name of output file")
output.add_argument("-m", "--multi_out", action="store_true", help="produces multiple output files")
# if -m specified, output file count = number of sequences in fasta file

args = parser.parse_args()
infile = args.infile
target_dir = args.target_directory
outfile = args.outfile

uids, coding_seqs = parse_fasta(infile)

# verified against ncbi 08Apr2019, plus Chris's exceptions in species.py
# allow for selenocysteine (TGA=U) (https://en.wikipedia.org/wiki/Selenocysteine)
# allow for pyrrolysine (TAG=O) (https://en.wikipedia.org/wiki/Pyrrolysine)
tt_11 = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': 'O',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'U', 'TGG': 'W',
}

# get list of dicts that count codon occurrences in cds
codon_dists = []
if args.multi_out:
    print(len(coding_seqs), "output files will be created")
    codon_counts = cd.get_protein_counts(coding_seqs, tt_11)

    for count_dict in codon_counts:
        # change counts dict into a frequency distribution
        codon_dist = cd.frequentize_counts(count_dict, tt_11)
        codon_dists.append(codon_dist)
        # print(len(codon_dists), codon_dists[0])
else:
    codon_counts = cd.get_org_counts(coding_seqs, tt_11)
    codon_dists.append(cd.frequentize_counts(codon_counts, tt_11))

try:
    os.mkdir(target_dir)
except OSError:
    pass

# uid and distribution indexing should still line up from original fasta parsing
# create csv file(s) to store distributions
for i, dist in enumerate(codon_dists):
    if args.multi_out:
        filename = os.path.join(target_dir, uids[i]) + '.csv'
    else:
        filename = os.path.join(target_dir, outfile)
    with open(filename, 'w') as fh:
        for codon, fraction in dist.items():
            fh.write(codon + ', ' + str(fraction) + '\n')



