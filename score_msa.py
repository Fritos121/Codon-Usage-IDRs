from Bio import AlignIO
from Bio.Alphabet import IUPAC, AlphabetEncoder

align_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_msa_codon3.txt'

msa_alphabet = AlphabetEncoder(IUPAC.ExtendedIUPACProtein(), '-.')

alignments = AlignIO.read(align_in, 'fasta', alphabet=msa_alphabet)
print(alignments)
quit()

# for i, alignment in enumerate(AlignIO.read(align_in, format, alphabet=msa_alphabet)):