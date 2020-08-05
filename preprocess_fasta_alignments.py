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
# print(msa_alphabet.letters)

# the two fed files should be in the same order (expected if using entire pipeline), allowing for easy indexing
coding_seqs = list(SeqIO.parse(cds_in, format, alphabet=IUPAC.IUPACUnambiguousDNA()))

codon_alignments = []
for i, alignment in enumerate(AlignIO.read(align_in, format, alphabet=msa_alphabet)):
    # don't need Seq object, just make string to make iterator more clear
    codon_seq = list(codon_iter(str(coding_seqs[i].seq)))
    codon_position = 0      # marks where we are in the list of codons
    translated_seq = ''     # translated from aa to DNA

    for aa in alignment.seq:
        if aa == '.' or aa == '-':
            translated_seq += aa*3

        else:
            codon = codon_seq[codon_position]
            # make sure codon isn't an open frame (at end of seq)
            if len(codon) % 3 == 0:
                # check if aa in aligned seq translates to corresponding codon in codon_seq
                trans_aa = try_translate(codon, tt_11)
                if aa.upper() == trans_aa and trans_aa is not None:
                    # check if aa in lower case
                    # if aa.islower():
                        # translated_seq += '...'   # . or -?
                    # else:
                    translated_seq += codon
                    codon_position += 1

                # else for when mismatch occurs?

            else:
                codon_position += 1

    # store SeqRecords later?
    #print(translated_seq)
    codon_alignments.append(translated_seq)
    # check that every codon was properly iterated over; if False, aa.upper() != output (via KeyError or legit mismatch)
    print(len(codon_seq) - 1 == codon_position)
    #break

print(len(codon_alignments))






