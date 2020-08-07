from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC, AlphabetEncoder
from codon_dist import codon_iter

def try_translate(codon, translation_table):
    try:
        aa = translation_table[codon]
    except KeyError:
        aa = None

    return aa


align_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_msa.txt'
cds_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_cds.fna'
format = 'fasta'

# no dots in file path until file extension
align_out = align_in.split('.', 1)[0] + '_codon.txt'

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

# check if alphabet required in long run... later
# adds gap both gap characters seen in PANTHER msa
msa_alphabet = AlphabetEncoder(IUPAC.ExtendedIUPACProtein(), '-.')

# the two fed files should be in the same order (expected if using entire pipeline), allowing for easy indexing
coding_seqs = list(SeqIO.parse(cds_in, format, alphabet=IUPAC.IUPACUnambiguousDNA()))

out_fh = open(align_out, 'w')
for i, alignment in enumerate(AlignIO.read(align_in, format, alphabet=msa_alphabet)):
    # coding seq id needs split, can be made more general for all fasta later...
    assert (alignment.id == coding_seqs[i].id.split(';')[0][4:]), "Alignment ID doesn't match coding sequence ID"
    #print(len(alignment.seq))
    #print(len(str(coding_seqs[i].seq)))

    # don't need Seq object, just make string to make iterator more clear
    codon_seq = list(codon_iter(str(coding_seqs[i].seq)))
    print(len(codon_seq))
    print(codon_seq)
    codon_position = 0      # marks where we are in the list of codons
    translated_seq = ''     # translated from aa to DNA
    j = 1   # counts number of non-gap aa; should end at len(codon_seq)... True will print if it does
    for aa in alignment.seq:

        if aa == '.' or aa == '-':
            translated_seq += aa * 3

        else:
            j += 1
            print(j, codon_position)
            codon = codon_seq[codon_position]
            #print(codon, codon_position * 3)
            # make sure codon isn't an open frame (at end of seq)
            if len(codon) % 3 == 0:
                # check if aa in aligned seq translates to corresponding codon in codon_seq
                trans_aa = try_translate(codon, tt_11)
                if aa.upper() == trans_aa and trans_aa is not None:
                    # check if aa in lower case
                    if aa.islower():
                        translated_seq += '...'
                    else:
                        translated_seq += codon

                    codon_position += 1

                # for when mismatch occurs; get alignment id  and alignment and cds positions
                else:
                    translated_seq += 'NNN'
                    codon_position += 1
                    #print('oops')

            else:
                codon_position += 1


    # check that every codon was properly iterated over; if False, aa.upper() != output (via KeyError or legit mismatch)
    print(j == len(codon_seq))
    #print(len(codon_seq) - 1 == codon_position)

    out_fh.write('>' + alignment.id + '\n')
    out_fh.write(translated_seq + '\n')

out_fh.close()






