from get_species_info_by_uid import get_protein_info
from preprocess_fasta_alignments import try_translate
from codon_dist import codon_iter
from vsl2 import run_vsl2b             # add to github?
import requests

# remove non-standard AAs so that vsl2 can run properly... use this tt for all in pipeline?
# TAA TAG and TGA all changed to A (which aa would be best?)
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

infile = r"D:\Orthologs\ecoli_gene_family_map_allortho.csv"

in_fh = open(infile, 'r')
uids = [x.split(',')[0] for x in in_fh.readlines()]
# uids = ['P0AAJ3', 'A0NAQ1']
in_fh.close()

# get list of uids
embl_accessions, tax_ids, ortho_mapping = get_protein_info(uids)

# embl_acc, uids, and tax_ids should have same length
# get CDS for gene
embl_base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
for uid, embl_acc in zip(uids, embl_accessions):
    if embl_acc is None:
        continue

    else:
        r = requests.get(embl_base_url + embl_acc)
        cds = ''.join(r.content.decode('utf-8').split("\n")[1:])
        # account for any http errors... will remove ortho from data for now; write uid info to std error?
        if 'status=' in cds or len(cds) == 0:
            print(uid, "removed from gene list due to error.")
            del ortho_mapping[uid]
        else:
            ortho_mapping[uid].append(cds)


disorder_scores = []
# use vsl to get disorder for each gene
for uid, info in ortho_mapping.items():
    # do not score any gene without a uid or cds
    # print(uid, info[2])
    if uid is None or info[2] is None:
        continue

    else:
        aa_seq = ''
        bad_aa = 0
        for i, codon in enumerate(codon_iter(info[2])):
            if len(codon) % 3 == 0:
                aa = try_translate(codon, tt_11)
                # dont add to seq fed into vsl2, but count
                if aa is None:
                    bad_aa += 1
                else:
                    aa_seq += aa

        disorder_list = run_vsl2b(aa_seq)[-1]

        # add count of bad aa so that proper fraction can be calculated
        fraction_disorder = disorder_list.count('D')/(len(disorder_list) + bad_aa)
        disorder_scores.append((uid, fraction_disorder))

# disorder_scores.sort(key=lambda x: x[1], reverse=True)    # no need to sort.. can sort in excel

with open(r'D:\Orthologs\ecoli_gene_disorder_ranked.csv', 'w') as fh:
    for uid, frac in disorder_scores:
        fh.write(uid + ',' + str(frac) + '\n')
