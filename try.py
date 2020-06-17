from Bio import Entrez
import re
from bioservices import UniProt

def get_match(pattern, search_str, index):
    match = re.search(pattern, search_str)
    try:
        if len(match.groups()) == 1:
            wanted_item = match.group(1)
        else:
            wanted_item = match.groups()

    except AttributeError:
        if pattern is refseq_pattern:
            missing_refseq.append(uids[index])
        elif pattern is taxid_pattern:
            missing_taxid.append(uids[index])
        elif pattern is ftp_pattern:
            missing_ftp.append([uids[index]])

    else:
        print(wanted_item)
        return wanted_item


def get_assembly_id(tx_id):

    search_term = 'txid{}'.format(tx_id)  # search term format used by NCBI Genome db
    handle1 = Entrez.esearch(db='genome', term=search_term) # gets internal ID for genome in genome db
    record1 = Entrez.read(handle1)
    #print(record1)
    org_id = record1['IdList'][0]    # what if there are more than 1... not encountered yet so don't worry
    #print(org_id)
    handle2 = Entrez.esummary(db='genome', id=org_id) # gives assembly of genome org
    record2 = Entrez.read(handle2)
    print(record2)
    assem_id = record2[0]["AssemblyID"]
    print(assem_id)

    # depending on design of program later, make this get entire record for assembly?

    return assem_id


Entrez.email = 'mlowry2@mymail.vcu.edu'
refseq_pattern = re.compile(r"DR\s+RefSeq.*(XM_.*).")
taxid_pattern = re.compile(r"OX\s+NCBI_TaxID=(\d+)")
ftp_pattern = re.compile(r"<FtpPath_RefSeq>(\S*(GCF_\S*))<")
missing_refseq = []
missing_taxid = []
missing_ftp = []

u = UniProt()
uids = ['A0NAQ1', 'A8DWE3'] # this will be gotten from bioservices.py somehow
prot_info = u.retrieve(uids, frmt='txt')

# get assembly ID

for i, record in enumerate(prot_info):
    if i == 0:
        refseq_ver = get_match(refseq_pattern, record, i)
        tax_id = get_match(taxid_pattern, record, i)

        assembly_id = get_assembly_id(tax_id)
        handle2 = Entrez.esummary(db='assembly', id=assembly_id, report='full')
        string = handle2.read()

        # gets ftp link for all cds in organism (download into <path>/tax_id/)
        ftp_matches = get_match(ftp_pattern, string, i)
        ftp_link = '/'.join(ftp_matches) + "_cds_from_genomic.fna.gz"
        print(ftp_link)







