# Get GC from Codon distribution
from Bio.SeqUtils import GC
from get_statsV3 import parse_fasta, parse_data
from statistics import mean, stdev, median


def codon_iter(seq):
    # takes DNA seq and yields generator of codons in seq
    # learn about generators and iterators
    for y in range(0, len(seq), 3):
        yield seq[y:y + 3]

def flip_trans_table(translation_table):
    # flip translation table so key:value is now AA:codon(s)
    tt_flip = {aa: [c for c, a in translation_table.items() if a == aa] for aa in set(translation_table.values())}

    return tt_flip

# make encapsulated function?
def get_Org_counts(protein_list, translation_table):
    '''
    gets codon and aa counts for whole organism
    param protein_list: list of each protein coding gene sequence in organism
    param translation_table: dictionary of 'codon':'AA' key-value pairs
    '''
    codon_counts = {}
    aa_counts = {}
    for x in protein_list:
        for codon in codon_iter(x):
            # don't include open frame end of gene
            if len(codon) % 3 == 0:
                # edit to check for ambigus bases with key 'ambig_base'
                # 4 total in .data file, skipped. make sure to catch when calc right when doing length based calc
                try:
                    aa = translation_table[codon]
                except KeyError:
                    continue

                try:
                    codon_counts[codon] += 1
                except KeyError:
                    codon_counts[codon] = 1

                try:
                    aa_counts[aa] += 1
                except KeyError:
                    aa_counts[aa] = 1

            else:
                print('Open frame')

    return codon_counts, aa_counts


def get_Prot_counts(protein_list, translation_table):
    '''
    makes list of dicts of aa_counts in protein from protein_list
    '''
    protein_counts = []
    for x in protein_list:
        aa_counts = {}
        for codon in codon_iter(x):
            # don't include open frame end of gene
            if len(codon) % 3 == 0:
                # remove
                try:
                    aa = translation_table[codon]
                except KeyError:
                    continue

                try:
                    aa_counts[aa] += 1
                except KeyError:
                    aa_counts[aa] = 1

        protein_counts.append(aa_counts)

    return protein_counts


def get_GC(codon_frac, aa_counts, translation_table):
    '''
    calculates GC for aa_counts dict provided (so for either org, or list of proteins)
    param codon_frac: dict of codon:prop of codon count to overall codon count for amino acid
    param aa_counts: dictionary of counts of each amino acid
    param translation_table: dictionary of 'codon':'AA' key-value pairs
    '''

    tt_flip = flip_trans_table(translation_table)
    gc_contrib = {}
    for aa, codons in tt_flip.items():
        if aa not in aa_counts:
            continue
        codon_gc = []
        for codon in codons:
            # not every codon in every protein, so some error handling... ew I know
            codon_gc.append((GC(codon) / 100) * (codon_frac[codon]))

        gc_contrib[aa] = sum(codon_gc) * (aa_counts[aa] / sum(aa_counts.values()))

    wanted_gc = (sum(gc_contrib.values()) * 100)

    return wanted_gc


def change_counts(codon_counts, translation_table):
    '''
    makes dict of codon:fraction this codon contributes to overall codon count for its amino acid
    param codon_counts: dictionary of codon counts in an organism
    '''

    tt_flip = flip_trans_table(translation_table)
    codon_frac = {}
    for aa, codons in tt_flip.items():
        # get total number of codons that translate to given aa
        total_codon = sum([codon_counts[codon] for codon in codons])
        for codon in codons:
            codon_frac[codon] = codon_counts[codon] / total_codon

    return codon_frac


if __name__ == '__main__':

    # verified against ncbi 08Apr2019, plus Chris's exceptins in species.py
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

    full_cds, cds_list, cds_gc = parse_fasta('CDS.fasta')

    codon_dict, aa_dict = get_Org_counts(cds_list, tt_11)
    frac_codon = change_counts(codon_dict, tt_11)

    protein_counts = get_Prot_counts(cds_list, tt_11)

    protein_gc = [get_GC(frac_codon, aac, tt_11) for aac in protein_counts]

    # print(protein_gc[0:2])

    print(GC(full_cds))
    print(mean(protein_gc))
    print(get_GC(frac_codon, aa_dict, tt_11))