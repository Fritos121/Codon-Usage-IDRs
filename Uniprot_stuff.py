from bioservices import UniProt
from Bio import Entrez
import re

def get_match(pattern, search_str, index):
    """
    Return the match(s) for either RegEx with error handling
    :param pattern: compiled RegEx obj
    :param search_str: string to search
    :param index: int to show current spot in list of proteins, for error handling only
    :return: either the wanted pattern found OR appends to list if error
    """

    match = re.search(pattern, search_str)
    try:
        wanted_item = match.group(1)
        print(wanted_item)
        return wanted_item

    except AttributeError:
        if pattern is refseq_pattern:
            missing_refseq.append(uids[index])
        elif pattern is taxid_pattern:
            missing_taxid.append(uids[index])

# maybe helper function in function that retrieves the assembly seqs to make genome
# not implemented in script 06Jun2020
def get_assembly_id(tx_id):
    """
    Used to get assembly ids for Chris code
    :param tx_id: string, NCBI tax_id
    :return: assembly ID for org
    """

    search_term = 'txid{}'.format(tx_id)  # search term format used by NCBI Genome db
    # clunky b/c of 2 api requests? only way I found atm
    handle1 = Entrez.esearch(db='genome', term=search_term)  # search term format used by NCBI Genome db
    record1 = Entrez.read(handle1)
    org_id = record1['IdList'][0]    # what if there are more than 1... not encountered yet so don't worry
    #print(org_id)
    handle2 = Entrez.esummary(db='genome', id=org_id)
    record2 = Entrez.read(handle2)
    #print(record2)
    assem_id = record2[0]['AssemblyID']

    return assem_id


Entrez.email = 'mlowry2@mymail.vcu.edu'
refseq_pattern = re.compile(r"DR\s+RefSeq.*(XM_.*).")
taxid_pattern = re.compile(r"OX\s+NCBI_TaxID=(\d+)")
missing_refseq = []
missing_taxid = []


# list of orthologus proteins
u = UniProt()
uids = ['A0NAQ1', 'A8DWE3'] # this will be gotten from bioservices.py somehow
prot_info = u.retrieve(uids, frmt='txt')

#print(prot_info[0], type(prot_info))
# loop through prot_info twice? once for protein's org info and another for protein info



# get protein coding DNA seq
for i, record in enumerate(prot_info):
    refseq_ver = get_match(refseq_pattern, record, i)
    tax_id = get_match(taxid_pattern, record, i)

    # get CDS for gene
    handle = Entrez.efetch(db='nuccore', id=refseq_ver, rettype='fasta', retmode='text')
    #if i == 1:
        #print(handle.readlines())
    CDS = ''.join([x.strip() for x in handle.readlines()[1:]])
    print(CDS)
    handle.close()

    # get CDSs of genome of org this gene found in and write codon dist to file in folder named 'tax_id'

    # get codon dist and compare to org codon dist in taxid folder?

if missing_refseq:
    print('Proteins Missing RefSeq version numbers: ' + ', '.join(missing_refseq))
if missing_taxid:
    print('Proteins Missing NCBI TaxIDs: ' + ', '.join(missing_taxid))






