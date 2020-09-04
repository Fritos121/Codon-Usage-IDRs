from Bio import AlignIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Alphabet import IUPAC, AlphabetEncoder
from codon_dist import flip_trans_table, codon_iter
from vsl2 import run_vsl2b
import re
import os


def fetch_org_distribution(taxonomy_id, base_dir):
    """
    Creates codon distribution for a given organism found in the given directory.
    Distributions files are expected to be in csv format, each line formatted as: codon, frequency
    :param taxonomy_id: string; NCBI tax_id
    :param base_dir: directory where folders named by tax_id are stored
    :return: dictionary of codons paired with their frequency in their aa codon distribution
    """
    import glob

    # add empty string to trick join to add one more set of filepath delimiters
    org_dir = os.path.join(base_dir, taxonomy_id, '')
    # grabs 1st csv file (should be only one if using whole pipeline)
    file_list = glob.glob(org_dir + '*_dist.csv')
    if len(file_list) != 1:
        exit('Too many or no csv files in directory. Ensure exactly one distribution file exists.')

    dist_fh = open(file_list[0], 'r')
    codon_dist = {codon: float(freq) for codon, freq in [line.strip().split(',') for line in dist_fh.readlines()]}

    return codon_dist


align_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_msa_codon3.txt'
source_org_dir = r'D:\Orthologs\Source_Org_Codon_Dist'
outdir = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792'

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

tt_flip = flip_trans_table(tt_11)

matrix = matlist.blosum62
# print(matrix)

msa_alphabet = AlphabetEncoder(IUPAC.ExtendedIUPACProtein(), '-.')
alignments = AlignIO.read(align_in, 'fasta', alphabet=msa_alphabet)

uid_pattern = re.compile(r'uid=(\S+?);')
tax_id_pattern = re.compile(r'tax_id=(\d+)')

# name file with uid of gene of interest as has been convention; expects it to be the first gene in the msa
uid = re.search(uid_pattern, alignments[0].id).group(1)

bad_codons = ['...', '---', 'NNN']    # codons that will not receive scores

# get disorder for each sequence in the msa
# make this a function later
disorder_info = []      # list of disorder score strings for each seq in msa (with gaps)
for alignment in alignments:
    aa_seq = ''
    bad_indexes = []
    for i, codon in enumerate(codon_iter(alignment.seq)):
        if codon in bad_codons:     # save location of bad codon
            bad_indexes.append(i)
        else:
            aa = tt_11[codon]
            aa_seq += aa

    disorder_list = list(run_vsl2b(aa_seq)[-1])     # last item returned is tuple of disorder score by residue

    # at these indexes, need to insert char to keep disorder scores registered with alignment
    for index in bad_indexes:
        disorder_list.insert(index, '-')       # insert - b/c . is used to indicate ordered residue

    disorder_scores = ''.join(disorder_list)
    disorder_info.append(disorder_scores)


out_fh = open(os.path.join(outdir, uid) + '_ortholog_msa_scores2.data', 'w')
out_fh.write("Identity,Percent Identity,Avg Blosum62 Score,Avg Frequency Score,Fraction Aligned,Fraction Disordered\n")

# calculate information for each column in alignment
for i in range(0, alignments.get_alignment_length(), 3):
    column = alignments[:, i:i+3]   # get MSA object with just the first column, size of one codon
    total_rows = len(column)
    aa_counts = {'x': 0}    # initialize with error aa

    running_blosum_score = 0
    running_freq_score = 0
    disorder_count = 0
    good_rows = 0       # count how many rows were used in blosum, freq, and disorder scoring

    for j, row in enumerate(column):
        codon1 = str(row.seq)
        aa1 = ''
        if codon1 in bad_codons:
            aa_counts['x'] += 1     # count error codon
            continue                # don't want to score bad codons

        else:
            aa1 = tt_11[codon1]
            try:
                aa_counts[aa1] += 1
            except KeyError:
                aa_counts[aa1] = 1

        good_rows += 1

        # get disorder score for codon; j = row, i = column in codon seq, i//3 = column in aa seq
        dis_score = disorder_info[j][i//3]
        if dis_score == 'D':
            disorder_count += 1

        # get frequency score for codon
        tax_id = re.search(tax_id_pattern, row.id).group(1)
        source_codon_dist = fetch_org_distribution(tax_id, source_org_dir)

        # 1 for frequent codons, 0 for rare/infrequent (when comparing the otho's source org codon dist to uniform dist)
        codon_count = len(tt_flip[aa1])
        if source_codon_dist[codon1] >= (1 / codon_count):
            running_freq_score += 1
        else:
            pass

        # get every row below current one in column
        for k in range(j + 1, total_rows):
            codon2 = str(column[k].seq)
            if codon2 in bad_codons:
                continue            # don't want to compare to bad codons; no score

            aa2 = tt_11[codon2]

            # all aa in matrix, but matrix isn't mirrored in dict, so try mirror of original key
            try:
                score = matrix[aa1, aa2]
            except KeyError:
                score = matrix[aa2, aa1]

            running_blosum_score += score

    # if no informational codons exist in column, or most common aa is an error
    identity = max(aa_counts, key=aa_counts.get)  # most common aa in column
    if good_rows == 0 or identity == 'x':
        out_fh.write("X,X,X,X,X,X\n")    # an X for every value recorded per column
        continue

    # calculate percent identity for column and fraction of column aligned properly
    percent_id = aa_counts[identity] / good_rows    # only want the identity of good columns
    fraction_aligned = good_rows / total_rows  # fraction of rows in column that fail to align

    # fraction of good rows that were disordered
    fraction_disordered = disorder_count / good_rows

    # calculate avg blosum62 and freq scores for column
    num_comparisons = (good_rows - 1) * good_rows / 2  # number of pairwise comparisons used to get column avg
    blosum_avg = running_blosum_score / num_comparisons
    avg_freq_score = running_freq_score / good_rows

    # use write_csv for clarity?
    out_fh.write(str(identity) + ',' + str(percent_id) + ',' + str(blosum_avg) + ',' + str(avg_freq_score) + ',' +
                 str(fraction_aligned) + ',' + str(fraction_disordered) + '\n')

out_fh.close()




