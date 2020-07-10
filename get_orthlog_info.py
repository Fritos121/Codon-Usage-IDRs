from bioservices import Panther
import re
import os

# where did i get allortholog file from panther (on ftp somewhere)?
# with a uniprot ID, read allortholog file to get panther family id that the protein is in

p = Panther()
# list of dicts; get persistant_ids for alignments we want to keep
family_msa = p.get_family_msa('PTHR42792')

# can get ortho info from Allorthologs file instead... didnt have it at time of making script
# interested in list of dicts using 'mapped' key; orthologs of given gene in given org
ortho = p.get_ortholog(['P04949'], '83333')

# if multiple genes from same org used, unmapped might be populated
if ortho['unmapped']:
    print(ortho['unmapped'])
else:
    print("All genes mapped.")

#print(family_msa[0], '\n')
#print(ortho['mapped'][0], '\n')

uniprotKB_id = re.compile(r".*?\|UniProtKB=(.*)")
ortho_mapping = {}  # maps ortholog persistent ids to its respective uid
print(len(family_msa))
for x in ortho['mapped']:
    match = re.search(uniprotKB_id, x['target_gene'])
    uid = match.group(1)
    pid = x['target_persistent_id']
    ortho_mapping[pid] = uid


# get alignments for all orthologs in family (all orthologs in family, but not all genes in family are orthologs)
orthos_msa = [alignment for alignment in family_msa if alignment['persistent_id'] in ortho_mapping.keys()]

print(len(orthos_msa), len(ortho_mapping))
# genes i've been working with ('A0NAQ1', 'A8DWE3')
# print(ortho_mapping['PTN002871129'], ortho_mapping['PTN001250064'])


# write alignment to file; placed in Panther family folder
base_dir = os.path.abspath(r"D:\Orthologs\Ortholog_Codon_Dist")

# if only one gene ever retrieved from family, take the gene name out of dir_name and add to file_name
dir_name = os.path.join(base_dir, "PTHR42792")
try:
    os.mkdir(dir_name)
except OSError:
    pass

file_name = os.path.join(dir_name, "P04949_ortholog_msa.txt")
with open(file_name, 'w') as fh:
    for alignment in orthos_msa:
        seq, pid = alignment.items()
        fh.write(">uid=" + ortho_mapping[pid[1]] + "\n")
        fh.write(seq[1] + "\n")

# pipe out uids for use in get_species_info_by_uid.py
# do this by making dict[<Ecoli_gene>] = [ortho_uids]?  then write out dict to file to save?
# dict[(<panther_family>, <Ecoli_gene>)] ???????
