from Bio import AlignIO, SeqIO

align_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_msa.txt'
cds_in = r'D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_cds.fna'
format = 'fasta'

alignments = AlignIO.parse(align_in, format)
coding_seqs = SeqIO.parse(cds_in, 'fasta')

for cds in coding_seqs:
    # test if seqs are same
    cds1 = cds.seq
    print(len(cds1))
    break

i=0
for alignment in alignments:
    print(type(alignment))
    print(alignment.get_alignment_length())
    print(alignment)
    #align1 = alignment.seq
    #print(len(align1))
    i += 1



