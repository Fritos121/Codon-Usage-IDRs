from get_species_info_by_uid import get_protein_info
from preprocess_fasta_alignments import try_translate
import requests
from Bio import Entrez

infile = r"D:\Orthologs\ecoli_gene_family_map_allortho.csv"

in_fh = open(infile, 'r')

uids = [x.split(',')[0] for x in in_fh.readlines()]

print(len(uids), uids)

# get list of uids
embl_accessions, tax_ids, ortho_mapping = get_protein_info(uids)

# embl_acc, uids, and tax_ids should have same length
# get CDS for gene
Entrez.email = "mlowry2@mymail.vcu.edu"
embl_base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
for uid, embl_acc in zip(uids, embl_accessions):
    if embl_acc is None:
        continue

    else:
        r = requests.get(embl_base_url + embl_acc)
        cds = ''.join(r.content.decode('utf-8').split("\n")[1:])
        # account for any http errors... will remove ortho from data for now; write uid info to std error?
        if 'status=' in cds:
            print('\n', uid, "removed from orthologs due to error.")
            del ortho_mapping[uid]
        else:
            ortho_mapping[uid].append(cds)


# use vsl to get disorder for each gene
for uid, info in ortho_mapping.items():
    # do not write any ortholog with incomplete information
    if uid is None or info[2] is None:
        continue

    else:
        print(info[2])
        quit()



