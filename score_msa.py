from Bio import AlignIO
from Bio.Alphabet import IUPAC, AlphabetEncoder
from codon_dist import codon_iter
import re
import glob
import os

align_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_msa_codon3.txt'
source_org_dir = r'D:\Orthologs\Source_Org_Codon_Dist'
outdir = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792'

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
    # add empty string to trick join to add one more set of filepath delimiters
    org_dir = os.path.join(source_org_dir, tax_id, '')
    # grabs 1st csv file (should be only one if using whole pipeline)
    source_dist_fh = open(glob.glob(org_dir + '*.csv')[0], 'r')
    source_codon_dist = {codon: float(freq) for codon, freq in [line.strip().split(',')
                                                                for line in source_dist_fh.readlines()]}

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

        # 1 if codon conserved from subject, 0 if codon was not
        if codon == subject_codons[i]:
            conserv_scores += '1'
        else:
            conserv_scores += '0'

        # 1 for frequent codons, 0 for rare/infrequent (when comparing the otho's source org codon dist to uniform dist)
        if source_codon_dist[codon] >= 0.5:
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



