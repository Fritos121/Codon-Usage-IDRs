from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC, AlphabetEncoder
from codon_dist import codon_iter

def try_translate(codon, translation_table):
    try:
        aa = translation_table[codon]
    except KeyError:
        aa = None

    return aa


if __name__ == '__main__':

    align_in = r'D:\Orthologs\Ortholog_Codon_Dist\2PTHR42792\P04949_ortholog_msa.txt'
    cds_in = r'D:\Orthologs\Ortholog_Codon_Dist\2PTHR42792\P04949_ortholog_cds.fna'
    error_in = r'D:\Orthologs\Ortholog_Codon_Dist\2PTHR42792\removed_uids.txt'
    format = 'fasta'

    align_out = ''.join(align_in.split('.')[:-1]) + '_codon.txt'

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

    # get list of uids that had errors; used to prevent assert errors
    error_fh = open(error_in, 'r')
    error_uids = [line.strip() for line in error_fh.readlines()]
    error_fh.close()

    # check if alphabet required in long run... later
    # adds gap both gap characters seen in PANTHER msa
    msa_alphabet = AlphabetEncoder(IUPAC.ExtendedIUPACProtein(), '-.')

    # the cds and alignment files should be in the same order except for uids in error file, allowing for easy indexing
    coding_seqs = list(SeqIO.parse(cds_in, format, alphabet=IUPAC.IUPACUnambiguousDNA()))

    codon_alignments = []
    found_flag = False
    i = 0
    for alignment in AlignIO.read(align_in, format, alphabet=msa_alphabet):
        # if uid is one that errored previously, then it was removed from the cds file
        if alignment.id in error_uids:
            continue  # don't increment i if error found; should realign files

        # coding seq header needs split, can be made more general for all fasta later...
        cds_header = coding_seqs[i].id
        cds_uid = cds_header.split(';')[0][4:]
        # any uids not in error list that don't match corresponding entry in cds file crash out
        assert (alignment.id == cds_uid), "Alignment ID: " + alignment.id + " doesn't match Coding Sequence ID: " + cds_uid

        #print(len(alignment.seq))
        #print(len(str(coding_seqs[i].seq)))

        # don't need Seq object, just make string to make iterator more clear
        codon_seq = list(codon_iter(str(coding_seqs[i].seq)))
        #print(len(codon_seq))
        #print(codon_seq)
        codon_position = 0      # marks where we are in the list of codons
        aligned_coding_seq = ''     # translated from aa to DNA
        for aa in alignment.seq:

            if aa == '.' or aa == '-':
                aligned_coding_seq += aa * 3

            else:
                #print(codon_position)
                # fix later; indicates an alignment seq that is longer than the coding seq
                try:
                    codon = codon_seq[codon_position]
                except IndexError:
                    print(alignment.id, cds_header)
                    aligned_coding_seq += 'QQQ'
                    continue      # counts how many AA longer align seq is than codon_iter
                    # break
                #print(codon, codon_position * 3)
                # make sure codon isn't an open frame (at end of seq)
                if len(codon) % 3 == 0:
                    # check if aa in aligned seq translates to corresponding codon in codon_seq
                    trans_aa = try_translate(codon, tt_11)
                    if aa.upper() == trans_aa and trans_aa is not None:
                        # check if aa in lower case
                        if aa.islower():
                            aligned_coding_seq += '...'
                        else:
                            aligned_coding_seq += codon

                    # for when mismatch occurs; get alignment id, and alignment and cds positions
                    else:
                        aligned_coding_seq += 'NNN'
                        # print('oops')

                codon_position += 1

        # check that every codon was properly iterated over; if False, aa.upper() != output (via KeyError or legit mismatch)
        # number codons (non-gap aa in align seq) = codon_pos + 1
        if not (codon_position == len(codon_seq) - 1):
            print(False)

        # make gene of interest that will have all genes aligned to it for analysis first
        if alignment.id in align_in:
            codon_alignments.insert(0, (cds_header, aligned_coding_seq))
            found_flag = True
        else:
            codon_alignments.append((cds_header, aligned_coding_seq))

        i += 1

    # don't write to file if not found
    if not found_flag:
        exit('Gene of interest not found.')

    with open(align_out, 'w') as fh:
        for header, seq in codon_alignments:
            fh.write('>' + header + '\n')
            fh.write(seq + '\n')




