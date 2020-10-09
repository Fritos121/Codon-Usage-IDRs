from Bio import AlignIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Alphabet import IUPAC, AlphabetEncoder
from codon_dist import flip_trans_table, codon_iter
import numpy as np
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


def calc_freq_score(aa, observed, flipped_tt):
    """
    Score a codon as either 1 or 0 for frequent and rare/infrequent, respectively. Uses uniform codon dist as expected.
    :param aa: amino codon translates to
    :param observed: frequency of codon from its source organism codon dist
    :return:
    """
    freq_observed = observed
    freq_expected = (1 / len(flipped_tt[aa]))  # uniform dist
    if freq_observed >= freq_expected:
        return 1
    else:
        return 0


if __name__ == '__main__':
    align_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR43553\P33941_ortholog_msa_codon.txt'
    source_org_dir = r'D:\Orthologs\Source_Org_Codon_Dist'
    outdir = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR43553'

    # don't want to remove non-standard aa here b/c tt_11 used for freq calc
    # making those codons A will change the distribution of Alanine performed later
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
    alignments = AlignIO.read(align_in, 'fasta', alphabet=msa_alphabet)     # need to make list when getting len()
    # print(len(alignments))   # number of orthologs use in final msa; another function in module will do this?

    uid_pattern = re.compile(r'uid=(\S+?);')
    tax_id_pattern = re.compile(r'tax_id=(\d+)')

    # name file with uid of gene of interest as has been convention; expects it to be the first gene in the msa
    uid = re.search(uid_pattern, alignments[0].id).group(1)

    bad_codons = ['...', '---', 'NNN']    # codons that will not receive scores

    # get disorder for each sequence in the msa
    # make this a function later
    disorder_strength = []        # list of numerical disorder strength score for each aa seq in msa (with gaps)
    disorder_letters = []       # list of disorder score letters for each seq in msa (with gaps)
    for alignment in alignments:
        aa_seq = ''
        bad_indexes = []
        for i, codon in enumerate(codon_iter(alignment.seq)):
            if codon in bad_codons:     # save location of bad codon
                bad_indexes.append(i)

            elif codon in ['TAA', 'TAG', 'TGA']:    # remove non-standard AAs so that vsl2 can run properly, make it A
                aa_seq += 'A'

            else:
                aa = tt_11[codon]
                aa_seq += aa

        # print(len(aa_seq), alignment.id)
        scores, letters = run_vsl2b(aa_seq)[2:4]    # last item returned is tuple of disorder score by residue
        letters = list(letters)      # needs to be list not tuple

        # at these indexes, need to insert char to keep disorder scores registered with alignment
        for index in bad_indexes:
            scores.insert(index, '-')       # insert - b/c . is used to indicate ordered residue
            letters.insert(index, '-')

        disorder_strength.append(scores)
        disorder_letters.append(letters)

    out_fh = open(os.path.join(outdir, uid) + '_ortholog_msa_scores3.data', 'w')
    out_fh.write("Identity,Percent Identity,Avg Blosum62 Score,Avg Frequency Score,Fraction Aligned,"
                 "Fraction Disordered,Avg Disorder Strength,Avg Frequency Ratio,Avg Log Odds Frequency\n")

    # calculate information for each column in alignment
    for i in range(0, alignments.get_alignment_length(), 3):
        column = alignments[:, i:i+3]   # get MSA object with just the first column, size of one codon
        total_rows = len(column)
        aa_counts = {'x': 0}    # initialize with error aa

        running_blosum_score = 0
        running_freq_score = 0
        freq_ratio_sum = 0
        freq_log_odds_sum = 0
        disorder_count = 0
        disorder_strength_sum = 0
        good_rows = 0       # count how many rows were used in scoring

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
            dis_letter = disorder_letters[j][i//3]
            if dis_letter == 'D':
                disorder_count += 1

            # same deal here but it is a float, add together to eventually get avg
            dis_score = disorder_strength[j][i//3]
            disorder_strength_sum += dis_score

            # get frequency score for codon
            tax_id = re.search(tax_id_pattern, row.id).group(1)
            source_codon_dist = fetch_org_distribution(tax_id, source_org_dir)

            # new freq score with expected for entire column given aa and codon dist
            observed_freq = source_codon_dist[codon1]
            expected_freq = sum([calc_freq_score(aa1, source_codon_dist[codon], tt_flip) * source_codon_dist[codon]
                                 for codon in tt_flip[aa1]])

            # 1 for frequent codons, 0 for rare/infrequent
            running_freq_score += calc_freq_score(aa1, observed_freq, tt_flip)

            # using observed / expected to as another score
            freq_ratio_sum += observed_freq / expected_freq

            # use log (ln) odds to get freq score
            freq_log_odds_sum += np.log(observed_freq / expected_freq)

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
            out_fh.write("X,X,X,X,X,X,X,X,X\n")    # an X for every value recorded per column
            continue

        # calculate percent identity for column and fraction of column aligned properly
        percent_id = aa_counts[identity] / good_rows    # only want the identity of good columns
        fraction_aligned = good_rows / total_rows  # fraction of rows in column that fail to align

        # calculate fraction of good rows that were disordered, and avg disorder strength of good rows in column
        fraction_disordered = disorder_count / good_rows
        avg_disorder_strength = disorder_strength_sum / good_rows

        # calculate avg blosum62 and freq scores for column
        num_comparisons = (good_rows - 1) * good_rows / 2  # number of pairwise comparisons used to get column avg
        blosum_avg = running_blosum_score / num_comparisons
        avg_freq_score = running_freq_score / good_rows
        avg_freq_ratio = freq_ratio_sum / good_rows
        avg_log_odds_freq = freq_log_odds_sum / good_rows

        # use write_csv for clarity?
        out_fh.write(str(identity) + ',' + str(percent_id) + ',' + str(blosum_avg) + ',' + str(avg_freq_score) + ',' +
                     str(fraction_aligned) + ',' + str(fraction_disordered) + ',' + str(avg_disorder_strength) + ',' +
                     str(avg_freq_ratio) + ',' + str(avg_log_odds_freq) + '\n')

    out_fh.close()




