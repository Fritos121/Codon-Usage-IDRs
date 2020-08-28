from Bio import AlignIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Alphabet import IUPAC, AlphabetEncoder
from codon_dist import flip_trans_table
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

# name file with uid of gene of interest as has been convention
uid = re.search(uid_pattern, alignments[0].id).group(1)

bad_codons = ['...', '---', 'NNN']    # codons that will not receive scores

# identity scores for columns
out_fh = open(os.path.join(outdir, uid) + '_ortholog_msa_ID_scores.data', 'w')
for i in range(0, alignments.get_alignment_length(), 3):
    column = alignments[:, i:i+3]   # get MSA object with just the first column, size of one codon
    total_rows = len(alignments)
    aa_counts = {'x': 0}    # initialize with error aa

    for row in column:
        # count error codon
        if row.seq in bad_codons:
            aa_counts['x'] += 1
            continue

        aa = tt_11[str(row.seq)]
        try:
            aa_counts[aa] += 1
        except KeyError:
            aa_counts[aa] = 1

    # print(aa_counts)
    keymax = max(aa_counts, key=aa_counts.get)  # most common aa in column
    identity = (keymax, aa_counts[keymax] / total_rows)

    # check if column with error as max is also the only codon at that position, deal with later
    if identity[0] == 'x' and identity[1] != 1.0:
        print(False)

    out_fh.write(str(identity) + ', ')

out_fh.close()


# use blosum62 matrix to score columns
out_fh = open(os.path.join(outdir, uid) + '_ortholog_msa_column_scores.data', 'w')
for i in range(0, alignments.get_alignment_length(), 3):
    column = alignments[:, i:i+3]
    column = [codon for codon in column if codon.seq not in bad_codons]     # filter out rows rows with no info
    total_rows = len(column)    # number of seqs in column

    # if no informational codons exist in column
    if total_rows == 0:
        column_info = ('X', 'X')
        out_fh.write(str(column_info) + ', ')
        continue

    num_comparisons = (total_rows - 1) * total_rows / 2  # number of pairwise comparisons used to get column avg
    running_score = 0
    running_freq_score = 0

    for j, row in enumerate(column):
        codon1 = str(row.seq)
        aa1 = tt_11[codon1]

        # get frequency score for codon
        header = row.id
        tax_id = re.search(tax_id_pattern, header).group(1)
        source_codon_dist = fetch_org_distribution(tax_id, source_org_dir)

        # 1 for frequent codons, 0 for rare/infrequent (when comparing the otho's source org codon dist to uniform dist)
        codon_count = len(tt_flip[aa1])
        if source_codon_dist[codon1] >= (1 / codon_count):
            running_freq_score += 1
        else:
            pass

        # get every alignment below current one in column
        for k in range(j+1, total_rows):
            codon2 = str(column[k].seq)
            aa2 = tt_11[codon2]

            # all aa in matrix, but matrix isn't mirrored in dict, so try mirror of original key
            try:
                score = matrix[aa1, aa2]
            except KeyError:
                score = matrix[aa2, aa1]
                
            running_score += score

    blosum_avg = running_score / num_comparisons
    avg_freq_score = running_freq_score / total_rows

    column_info = (blosum_avg, avg_freq_score)
    out_fh.write(str(column_info) + ', ')

out_fh.close()


