from Bio import AlignIO
from Bio.Alphabet import IUPAC, AlphabetEncoder
from codon_dist import codon_iter, flip_trans_table
import re
import os

def fetch_org_distribution(taxonomy_id, base_dir):
    """
    Creates codon distribution for a given organism found in the given directory
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

msa_alphabet = AlphabetEncoder(IUPAC.ExtendedIUPACProtein(), '-.')

alignments = list(AlignIO.read(align_in, 'fasta', alphabet=msa_alphabet))
# gene of interest expected to be first gene in msa file (should be if preprocessed)
subject = alignments.pop(0)
subject_codons = list(codon_iter(str(subject.seq)))

uid_pattern = re.compile(r'uid=(\S+?);')
tax_id_pattern = re.compile(r'tax_id=(\d+)')

# name file with uid of gene of interest as has been convention
uid = re.search(uid_pattern, subject.id).group(1)
out_fh = open(os.path.join(outdir, uid) + '_ortholog_msa_scores.data', 'w')

bad_codons = ['...', '---', 'NNN']    # codons that will not receive scores

for alignment in alignments:
    header = alignment.id
    # tax_id needed to obtain files containing source org codon dist... used to obtain frequency score for codon
    tax_id = re.search(tax_id_pattern, header).group(1)
    source_codon_dist = fetch_org_distribution(tax_id, source_org_dir)

    alignment_codons = list(codon_iter(str(alignment.seq)))
    freq_scores = ''
    conserv_scores = ''

    # alignments all same length, so enumerate can be used
    for i, codon in enumerate(alignment_codons):
        # placeholder no-score value for gaps and ambiguous codons
        #if codon in bad_codons or subject_codons[i] in bad_codons:
        if codon in bad_codons:     # assumes both msa's gaps are 1:1; tested, true so far
            conserv_scores += 'x'
            freq_scores += 'x'
            continue

        # 1 if aa is conserved from subject (synonymous codons), 0 if not
        aa = tt_11[codon]
        if aa == tt_11[subject_codons[i]]:
            conserv_scores += '1'
        else:
            conserv_scores += '0'

        # 1 for frequent codons, 0 for rare/infrequent (when comparing the otho's source org codon dist to uniform dist)
        codon_count = len(tt_flip[aa])      # counts number of codons per aa to get uniform dist below
        if source_codon_dist[codon] >= (1/codon_count):
            freq_scores += '1'
        else:
            freq_scores += '0'

    out_fh.write('>' + header + '\n')
    out_fh.write(str(alignment.seq) + '\n')
    out_fh.write(conserv_scores + '\n' + freq_scores + '\n')


    #print(conserv_scores)
    print(conserv_scores.count('1')/(conserv_scores.count('1') + conserv_scores.count('0')))
    #print(freq_scores)
    #break

out_fh.close()



