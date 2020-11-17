from get_species_info_by_uid import get_protein_info
from preprocess_fasta_alignments import try_translate
from codon_dist import codon_iter
from vsl2 import run_vsl2b
import requests
import sys


def parse_fasta(file_handle):
    """
    Separate fasta headers and sequences into two lists
    :param file_handle: file handle of fasta file
    :return: list of headers, list of sequences
    """
    entry_headers = []
    entry_seqs = []
    for line in file_handle:
        line = line.strip()
        if line.startswith('>'):
            entry_headers.append(line)

        else:
            sequence = line
            entry_seqs.append(sequence)

    return entry_headers, entry_seqs


def get_disorder_scores(coding_seq, trans_table):
    """
    Translates DNA sequence to proteins; uses VSL2 to predict disorder in protein sequence.
    :param coding_seq: string; DNA protein coding sequence
    :param trans_table: dictionary of 'codon':'AA' key-value pairs
    :return: disorder scores from VSL2, and the length of the analyzed sequence
    """
    aa_seq = ''
    for i, codon in enumerate(codon_iter(coding_seq)):
        if len(codon) % 3 == 0:
            aa = try_translate(codon, trans_table)
            # dont add to seq fed into vsl2
            if aa is None:
                continue
            else:
                aa_seq += aa

    # -1 gets hard score, -2 gets soft score
    scores = run_vsl2b(aa_seq)[-2]
    seq_length = len(aa_seq)

    return scores, seq_length


if __name__ == '__main__':

    # ensure proper command line arguments are passed.
    if len(sys.argv) != 3:
        exit(f"Required positional arguments: {sys.argv[0]} <input_file> <output_file>")

    infile = sys.argv[1]
    outfile = sys.argv[2]

    headers = []
    seqs = []

    # program can take either a fasta file or txt file
    if infile.endswith('.txt'):
        in_fh = open(infile, 'r')
        uids = [line.strip() for line in in_fh.readlines()]
        in_fh.close()
        ortho_mapping = get_protein_info(uids)

        # get coding seqs from ENA
        embl_base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
        bad_ids = []
        for header, info in ortho_mapping.items():
            embl_acc = info[0]
            # dont want any gene where a coding sequence cannot be retrieved
            if embl_acc is None:
                continue

            else:
                r = requests.get(embl_base_url + embl_acc)
                cds = ''.join(r.content.decode('utf-8').split("\n")[1:])
                # account for any http errors... will remove ortho from data for now; write uid info to std error?
                if 'status=' in cds or len(cds) == 0:
                    print(header, "removed from gene list due to error in retrieving coding sequence.")
                    continue

                else:
                    headers.append(header)
                    seqs.append(cds)

    elif infile.endswith('.fasta') or infile.endswith('.fna'):
        in_fh = open(infile, 'r')
        headers, seqs = parse_fasta(in_fh)
        in_fh.close()

    else:
        exit("Please provide a text file with UniProt IDs separated by newlines, or a DNA fasta file.")

    # remove non-standard AAs so that vsl2 can run properly.
    # TAA TAG and TGA all changed to A (neutral disorder profile)
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
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'A', 'TAG': 'A',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'A', 'TGG': 'W',
    }

    with open(outfile, 'w') as out_fh:
        for header, seq in zip(headers, seqs):
            dis_scores, aa_seq_len = get_disorder_scores(seq, tt_11)

            # divide by length of amino acid seq used in vsl2 to get proper avg
            # remember that failed translations are not scored in function
            avg_disorder_strength = sum(dis_scores) / aa_seq_len
            out_fh.write(header + ',' + str(avg_disorder_strength) + '\n')
