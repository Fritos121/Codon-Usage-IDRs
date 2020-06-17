from bioservices import Panther
#from bioservices import UniProt
import re


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
ortho_pids = []
ortho_uids = []
print(len(family_msa))
for x in ortho['mapped']:
    match = re.search(uniprotKB_id, x['target_gene'])
    uid = match.group(1)
    pid = x['target_persistent_id']
    ortho_pids.append(pid)
    ortho_uids.append(uid)


# get alignments for all orthologs in family (all orthologs in family, but not all genes in family are orthologs)
orthos_msa = [alignment for alignment in family_msa if alignment['persistent_id'] in ortho_pids]

# can remove if not needed, just helps
ortho_ids = list(zip(ortho_pids, ortho_uids))
print(ortho_ids[0:2])

print(len(orthos_msa), len(ortho_pids))

# now, with uniprot IDs of orthologs, we need to get to the mRNA or DNA seq at NCBI though uniprot

