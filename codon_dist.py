# Get GC from Codon distribution
from Bio.SeqUtils import GC
from get_stats import parse_fasta, parse_data
from statistics import mean


def codon_iter(seq):
    """
    Takes DNA seq and yields generator of codons in seq
    :param seq: string; DNA sequence
    """

    for y in range(0, len(seq), 3):
        yield seq[y:y + 3]


def flip_trans_table(translation_table):
    """
    Flip translation table so key:value is now AA:codon(s)
    :param translation_table: dictionary of 'codon':'AA' key-value pairs
    :return: dictionary of flipped translation table
    """

    tt_flip = {aa: [c for c, a in translation_table.items() if a == aa] for aa in set(translation_table.values())}

    return tt_flip


def get_org_counts(protein_list, translation_table):
    """
    Gets codon counts for whole organism/list of given proteins
    :param protein_list: list of each protein coding gene sequence in organism
    :param translation_table: dictionary of 'codon':'AA' key-value pairs
    :return dictionary of codon counts
    """

    codon_counts = {codon: 0 for codon in translation_table.keys()}
    # for any codon present in org but not present in the passed translation table
    codon_counts['NNN'] = 0
    for x in protein_list:
        for codon in codon_iter(x):
            # don't include open frame end of gene
            if len(codon) % 3 == 0:
                try:
                    codon_counts[codon] += 1
                except KeyError:
                    codon_counts['NNN'] += 1

    return codon_counts


def get_protein_counts(protein_list, translation_table):
    """
    Makes list of dicts of codon counts in protein from protein_list
    :param protein_list: list of protein coding DNA sequences
    :param translation_table: dictionary of 'codon':'AA' key-value pairs
    :return list of dictionaries, one dict for each protein in protein_list
    """

    protein_counts = []
    for x in protein_list:
        codon_counts = {codon: 0 for codon in translation_table.keys()}
        # for any codon present in org but not present in the passed translation table
        codon_counts['NNN'] = 0
        for codon in codon_iter(x):
            # don't include open frame end of gene
            if len(codon) % 3 == 0:
                try:
                    codon_counts[codon] += 1
                except KeyError:
                    codon_counts['NNN'] += 1

        protein_counts.append(codon_counts)

    return protein_counts


def translate_counts(codon_counts, translation_table):
    """
    Translates a list of codon count dictionaries using passed translation table
    :param codon_counts: list of dicts of codon counts
    :param translation_table: dictionary of 'codon':'AA' key-value pairs
    :return: a list of dicts of amino acid counts
    """

    translated_counts = []
    for count_dict in codon_counts:
        aa_counts = {aa: 0 for aa in set(translation_table.values())}
        # for any codon present in count dict but not present in the passed translation table
        aa_counts['x'] = 0

        for codon, count in count_dict.items():
            try:
                aa = translation_table[codon]
                aa_counts[aa] += count
            except KeyError:
                aa_counts['x'] += count

        translated_counts.append(aa_counts)

    return translated_counts


def get_GC(codon_frac, aa_counts, translation_table):
    """
    Calculates GC for aa_counts dict provided (so for either org, or list of proteins)
    :param codon_frac: dict of codon:prop of codon count to overall codon count for amino acid
    :param aa_counts: dictionary of counts of each amino acid
    :param translation_table: dictionary of 'codon':'AA' key-value pairs
    :return int; calculated GC
    """

    tt_flip = flip_trans_table(translation_table)
    gc_contrib = {}
    for aa, codons in tt_flip.items():
        if aa not in aa_counts:
            continue
        codon_gc = []
        for codon in codons:
            # not every codon in every protein, so some error handling
            codon_gc.append((GC(codon) / 100) * (codon_frac[codon]))

        gc_contrib[aa] = sum(codon_gc) * (aa_counts[aa] / sum(aa_counts.values()))

    wanted_gc = (sum(gc_contrib.values()) * 100)

    return wanted_gc


def frequentize_counts(codon_counts, translation_table):
    """
    Makes dict of codon:<fraction this codon contributes to overall codon count for its amino acid>
    :param codon_counts: dictionary of codon counts
    :param translation_table: dictionary of 'codon':'AA' key-value pairs
    :return dictionary of codon distributions per their respective amino acid
    """

    tt_flip = flip_trans_table(translation_table)
    codon_frac = {}
    for aa, codons in tt_flip.items():
        # get total number of codons that translate to given aa
        total_codon = sum([codon_counts[codon] for codon in codons])
        for codon in codons:
            # since all codons from trans_table initialized, some might add to zero total
            try:
                codon_frac[codon] = codon_counts[codon] / total_codon
            except ZeroDivisionError:
                codon_frac[codon] = 0

    return codon_frac


if __name__ == '__main__':

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

    full_cds, cds_list, cds_gc = parse_fasta('CDS.fasta')

    codon_dict = get_org_counts(cds_list, tt_11)
    aa_dict = translate_counts([codon_dict], tt_11).pop(0)
    frac_codon = frequentize_counts(codon_dict, tt_11)

    protein_counts = get_protein_counts(cds_list, tt_11)
    protein_aa = translate_counts(protein_counts, tt_11)

    protein_gc = [get_GC(frac_codon, aac, tt_11) for aac in protein_aa]

    # print(protein_gc[0:2])

    print(GC(full_cds))
    print(mean(protein_gc))
    print(get_GC(frac_codon, aa_dict, tt_11))