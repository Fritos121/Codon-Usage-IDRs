import sys
import os
import codon_dist as cd


def parse_fasta(filename):
    with open(filename, 'r') as fh:
        uids = []
        cds_list = []
        for line in fh:
            if line.startswith('>'):
                uid = line.split(';')[0][5:]
                uids.append(uid)

            else:
                seq = line.strip()
                cds_list.append(seq)

    return uids, cds_list


# fasta file created by get_species_info_from_uid.py
infile = sys.argv[1]

# going to be the family folder + /Codon_Distributions
# make default directory same dir as input file (build from it?), and make this optional arg?
target_dir = sys.argv[2]
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

codon_counts = cd.get_Prot_counts(coding_seqs, tt_11, mode='codon')

codon_dists = []
i=0
for count_dict in codon_counts:
    codon_dist = cd.change_counts(count_dict, tt_11)
    codon_dists.append(codon_dist)
    if i == 0:
        print(len(codon_dists), codon_dists[0])
        quit()

try:
    os.mkdir(target_dir)
except OSError:
    pass








