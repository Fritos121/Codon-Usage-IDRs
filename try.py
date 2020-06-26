from Bio import Entrez
from time import sleep
import re


def get_protein_info(uniprot_ids):
    from bioservices import UniProt

    refseq_vers = []
    taxonomy_ids = []
    missing_refseq = []
    missing_taxid = []
    u = UniProt()

    uniprot_records = u.retrieve(uniprot_ids, frmt='txt')
    # print(prot_info[0])

    for i, record in enumerate(uniprot_records):
        refseq_ver = get_match(refseq_pattern, record, missing_refseq, uniprot_ids, i)
        taxonomy_id = get_match(taxid_pattern, record, missing_taxid, uniprot_ids, i)

        # make dict to connect ortho uid with refseq ver?
        refseq_vers.append(refseq_ver)  # refseq version number for coding seq
        taxonomy_ids.append(taxonomy_id)  # tax_id of organism protein belongs to

    if missing_refseq:
        print('Proteins Missing RefSeq version numbers: ' + ', '.join(missing_refseq))
    if missing_taxid:
        print('Proteins Missing NCBI TaxIDs: ' + ', '.join(missing_taxid))

    return refseq_vers, taxonomy_ids


# keep generic function even though on len() == 1 used in main()
def get_match(pattern, search_str, missing_list, id_list, index):
    match = re.search(pattern, search_str)
    try:
        if len(match.groups()) == 1:
            wanted_item = match.group(1)
        else:
            wanted_item = match.groups()

    except AttributeError:
        missing_list.append(id_list[index])
        # add placeholder value to keep indexing between refseq and tax_id 1:1?
        # return None

    else:
        print(wanted_item)
        return wanted_item


# set email in main()
def get_assembly_id(tx_id):

    search_term = 'txid{}'.format(tx_id)  # search term format used by NCBI Genome db
    handle1 = Entrez.esearch(db='genome', term=search_term) # gets internal ID for genome in genome db
    record1 = Entrez.read(handle1)
    # print(record1)
    org_id = record1['IdList'][0]    # what if there are more than 1... not encountered yet so don't worry
    # print(org_id)
    handle2 = Entrez.esummary(db='genome', id=org_id) # gives assembly of genome org
    record2 = Entrez.read(handle2)
    # print(record2)
    assem_id = record2[0]["AssemblyID"]
    # print(assem_id)

    # depending on design of program later, make this get entire report for assembly?

    return assem_id


if __name__ == "__main__":
    Entrez.email = 'mlowry2@mymail.vcu.edu'
    refseq_pattern = re.compile(r"DR\s+RefSeq.*(XM_.*).")
    taxid_pattern = re.compile(r"OX\s+NCBI_TaxID=(\d+)")
    ftp_pattern = re.compile(r"<FtpPath_RefSeq>\S+(/genomes\S+GCF_\S+)<")

    uids = ['A0NAQ1', 'A8DWE3']  # this will be gotten from bioservices.py (piped from command line)
    refseq_versions, tax_ids = get_protein_info(uids)
    # print(refseq_versions, tax_ids)
    # if implemented, remove placeholders
    unique_tax_ids = list(set(tax_ids))

    # get partial ftp links to be used by ftp module
    ftp_partials = []  # list of partial ftp links to be used later
    missing_ftp = []
    for i, tax_id in enumerate(unique_tax_ids):
        assembly_id = get_assembly_id(tax_id)
        handle = Entrez.esummary(db='assembly', id=assembly_id, report='full')
        string = handle.read()
        # print(string)

        # gets ftp (partial) link for all cds in organism
        ftp_link = get_match(ftp_pattern, string, missing_ftp, unique_tax_ids, i)
        ftp_partials.append((ftp_link, tax_id))
        sleep(0.4)  # prevent timeouts from ncbi

    print(ftp_partials)
    if missing_ftp:
        print('Proteins Missing FTP Link for Assembly: ' + ', '.join(missing_ftp))

    # do something with ortho CDS here. Determine how to store
    # refseq_versions, uids, and tax_ids should have same length

    # for ortho in enumerate(refseq_versions):

    # pass match object into function...?
