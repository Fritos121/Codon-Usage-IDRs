from bioservices import Panther
import re
import os
import sys
import platform
import argparse

# E. coli genes need disorder, but need some ordered controls too
# make txt file for list of uids of the above genes, grep allortholog file for the families

# where did i get allortholog file from panther (on ftp somewhere)?
# with a uniprot ID, read allortholog file to get panther family id that the protein is in


# command line arguments
# make it take uniprot gene id, panther family id, and taxid first, then try to make it take csv file with all info

parser = argparse.ArgumentParser()
parser.add_argument("target_dir", help="directory to store output file(s). Outfile name produced automatically.")
#parser.add_argument("-o", "--outfile", type=argparse.FileType('w'), help="name of output MSA file")
gene_info = parser.add_mutually_exclusive_group(required=True)
gene_info.add_argument("-g", "--gene_info", nargs=3, help="One gene's UniProt ID, PANTHER Family ID, and Taxonomy ID, in that order")
gene_info.add_argument("-i", "--infile", type=argparse.FileType('r'),
                       help="CSV file with one or more sets of UniProt IDs, PANTHER Family IDs, and Taxonomy IDs, in that order")


args = parser.parse_args()
base_dir = args.target_dir
#out_name = args.outfile
gene_ids = args.gene_info
in_fh = args.infile

print(base_dir)
#print(out_name)
print(gene_ids)
print(in_fh)

if in_fh:
    genes = [line.strip().split(',') for line in in_fh]

# gene_info was passed if in_fh wasn't
else:
    genes = []
    genes.append(gene_ids)
    

p = Panther()
for gene_list in genes:
    gene = gene_list[0]
    family = gene_list[1]
    tax_id = gene_list[2]
    print(gene, family, tax_id)
    quit()

# list of dicts; get persistent_ids for alignments we want to keep
family = 'PTHR43553'
family_msa = p.get_family_msa(family)

# interested in list of dicts using 'mapped' key; orthologs of given gene in given org
gene = 'P33941' # ['P04949']
ortho = p.get_ortholog(gene, '83333')

# if multiple genes from same org used, unmapped might be populated
if ortho['unmapped']:
    print(ortho['unmapped'])
else:
    print("All genes mapped.")

# print(family_msa[0], '\n')
# print(ortho['mapped'][0], '\n')

uniprotKB_id = re.compile(r".*?\|UniProtKB=(.*)")
ortho_mapping = {}  # maps ortholog persistent ids to its respective uid
print(len(family_msa))

missing_target_persis_id = []
for x in ortho['mapped']:
    # print(x)
    match = re.search(uniprotKB_id, x['target_gene'])
    uid = match.group(1)

    # can be obtained by searching target_gene manually in Panther db... go to end of tree msa
    # make function if target_gene or id have same issue. make error out if exceeding certain number of errors?
    try:
        pid = x['target_persistent_id']
    except KeyError:
        # pid = input("\nSearch " + uid + " in Panther db and enter Persistent ID from end of Tree: ")
        missing_target_persis_id.append(uid)
        continue

    ortho_mapping[pid] = uid

# adds the gene of interest to ortho_mapping so its information is included in later analysis and processing
else:
    uid = x['id']
    pid = x['persistent_id']
    ortho_mapping[pid] = uid

if missing_target_persis_id:
    print('\n{} Ortholog(s) Removed Due To Missing Target Persistent ID: '.format(len(missing_target_persis_id))
          + ', '.join(missing_target_persis_id))

# print(ortho_mapping)

# get alignments for all orthologs in family (all orthologs in family, but not all genes in family are orthologs)
orthos_msa = [alignment for alignment in family_msa if alignment['persistent_id'] in ortho_mapping.keys()]

print(len(orthos_msa), len(ortho_mapping))
# genes i've been working with ('A0NAQ1', 'A8DWE3')
# print(ortho_mapping['PTN002871129'], ortho_mapping['PTN001250064'])


# write alignment to file; placed in Panther family folder; allows windows and linux file paths
# once testing completed, make program take sys arg for base_dir, remove this
if platform.system() == "Windows":
    base_dir = r"D:\Orthologs\Ortholog_Codon_Dist"
elif platform.system() == "Linux":
    base_dir = os.path.abspath("/mnt/d/Orthologs/Ortholog_Codon_Dist")

# base_dir = sys.argv[1]
fn_base = "_ortholog_msa.txt"

# if only one gene ever retrieved from family, take the gene name out of dir_name and add to file_name
dir_name = os.path.join(base_dir, family)
try:
    os.mkdir(dir_name)
except OSError:
    pass

file_name = os.path.join(dir_name, gene) + fn_base
with open(file_name, 'w') as fh:
    for alignment in orthos_msa:
        seq, pid = alignment.items()
        fh.write(">" + ortho_mapping[pid[1]] + "\n")
        fh.write(seq[1] + "\n")

