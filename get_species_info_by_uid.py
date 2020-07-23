# added to prevent error (comment out to see it)
# from gevent import monkey as curious_george
# curious_george.patch_all(thread=False, select=False)

from Bio import Entrez
from time import sleep
import re
import requests
import sys


def get_protein_info(uniprot_ids):
    """
    Retrieves NCBI RefSeq version numbers and taxonomy ids for list of proteins. Creates a dict to map each protein's
    uid to its refseq version and tax id
    :param uniprot_ids: List of Uniprot IDs, e.g., "A5HYD5"
    :return: List of Refseq Version numbers, list of taxonmy ids, and dictionary mapping each uid to its info
    """
    from bioservices import UniProt

    embl_accs = []
    taxonomy_ids = []
    # refseq_vers = []
    missing_embl = []
    missing_taxid = []
    # missing_refseq = []
    orthos_map = {}
    u = UniProt()

    # uniprot_records = list(map(lambda x: x.decode("utf-8"), u.retrieve(uniprot_ids, frmt='txt')))  # WSL cmdline
    uniprot_records = u.retrieve(uniprot_ids, frmt='txt')  # PyCharm

    embl_pattern = re.compile(r"DR\s+EMBL;.*?;\s+(.*?);")
    taxid_pattern = re.compile(r"OX\s+NCBI_TaxID=(\d+)")
    # refseq_pattern = re.compile(r"DR\s+RefSeq.*(XM_.*).")

    for i, record in enumerate(uniprot_records):
        # record = "fffffffffffffffffff"    # test placeholder Nones

        embl_acc = get_match(embl_pattern, record, uniprot_ids, missing_embl, i)
        taxonomy_id = get_match(taxid_pattern, record, uniprot_ids, missing_taxid, i)
        # refseq_ver = get_match(refseq_pattern, record, uniprot_ids, missing_refseq, i)

        embl_accs.append(embl_acc)  # EMBL accession number for coding seq
        taxonomy_ids.append(taxonomy_id)  # tax_id of organism protein belongs to
        # refseq_vers.append(refseq_ver)  # refseq version number for coding seq
        orthos_map[uniprot_ids[i]] = [embl_acc, taxonomy_id]  # map protein info to its uid

    if missing_embl:
        print('\n{} Proteins Missing EMBL Accession Numbers: '.format(len(missing_embl)) + ', '.join(missing_embl))

    if missing_taxid:
        print('\n{} Proteins Missing NCBI TaxIDs: '.format(len(missing_taxid)) + ', '.join(missing_taxid))

    # if missing_refseq:
    #   print('\n{} Proteins Missing RefSeq Version Numbers: '.format(len(missing_refseq)) + ', '.join(missing_refseq))

    return embl_accs, taxonomy_ids, orthos_map


# keep generic function even though only len() == 1 used in main()
def get_match(pattern, search_str, id_list, missing_list, index):
    """
    Find first match using RegEx; Allows for tracking of what IDs can't get the specific pattern match and placeholders
    to maintain indexing
    :param pattern: compiled RegEx obj
    :param search_str: string to search
    :param id_list: list of identifiers
    :param missing_list: list to add identifier if pattern not found
    :param index: int; to show current spot in list of ids, for error handling only
    :return: Either the wanted pattern found (string), OR appends to list if error and returns placeholder None
    """

    match = re.search(pattern, search_str)
    try:
        if len(match.groups()) == 1:
            wanted_item = match.group(1)
        else:
            wanted_item = match.groups()

    except AttributeError:
        missing_list.append(id_list[index])
        # add placeholder value to keep indexing between refseq and uid 1:1, zipped later
        return None

    else:
        # print(wanted_item)
        return wanted_item


# set email in main()
def get_assembly_id(tx_id):
    """
    Retrieves NCBI Assembly ID for a given NCBI Taxonomy ID
    :param tx_id: string; NCBI tax id
    :return: string; NCBI Assembly ID
    """
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

    return assem_id


if __name__ == "__main__":
    '''
    # ensure proper command line arguments are passed.
    if len(sys.argv) != 3:
        exit("Required positional arguments: {} <infile> <outfile>".format(sys.argv[0]))

    infile = sys.argv[1]
    outfile = sys.argv[2]

    in_fh = open(infile)

    # get list of uids
    uids = [x.strip().replace(">", "") for x in in_fh.readlines() if x.startswith(">")]
    '''
    Entrez.email = 'mlowry2@mymail.vcu.edu'

    uids = ['A0NAQ1', 'A8DWE3']
    embl_accessions, tax_ids, ortho_mapping = get_protein_info(uids)
    # refseq_versions, tax_ids, ortho_mapping = get_protein_info(uids)
    # print(refseq_versions, tax_ids)

    # remove placeholder Nones
    tax_ids = [x for x in tax_ids if x is not None]
    unique_tax_ids = list(set(tax_ids))

    # get partial ftp links to be used by get_assembly_seqs_from_ftp.py
    ftp_pattern = re.compile(r"<FtpPath_RefSeq>\S+(/genomes\S+GCF_\S+)<")
    ftp_partials = []
    missing_ftp = []
    for i, tax_id in enumerate(unique_tax_ids):
        assembly_id = get_assembly_id(tax_id)
        handle = Entrez.esummary(db='assembly', id=assembly_id, report='full')
        # string = handle.read().decode("utf-8")  # for WSL
        string = handle.read()
        # print(string)

        # gets ftp (partial) link for all cds in organism
        ftp_link = get_match(ftp_pattern, string, missing_ftp, unique_tax_ids, i)
        ftp_partials.append((ftp_link, tax_id))
        sleep(0.4)  # prevent timeouts from ncbi

    if missing_ftp:
        print('\n{} Proteins Missing FTP Link for Assembly: '.format(len(missing_ftp)) + ', '.join(missing_ftp))

    # refseq_versions or embl_acc, uids, and tax_ids should have same length
    # get CDS for gene
    embl_base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
    for uid, embl_acc in zip(uids, embl_accessions):
        if embl_acc is None:
            continue

        else:
            r = requests.get(embl_base_url + embl_acc)
            cds = ''.join(r.content.decode('utf-8').split("\n")[1:])

            # use if refseq used to aquire cds
            # handle = Entrez.efetch(db='nuccore', id=refseq_ver, rettype='fasta', retmode='text')
            # cds = ''.join([x.strip() for x in handle.readlines()[1:]])
            # handle.close()

            ortho_mapping[uid].append(cds)


    # write cds to file
    outfile = r"D:\Orthologs\Ortholog_Codon_Dist\PTHR42792\P04949_ortholog_cds2.fasta"
    with open(outfile, 'w') as fh:
        for uid, info in ortho_mapping.items():
            if uid is None or None in info:
                continue

            # use either ";refseq_ver=" or ";embl_acc=" for info[0] depending on what was used above
            fh.write(">uid=" + uid + ";embl_acc=" + info[0] + ";tax_id=" + info[1] + "\n")
            fh.write(info[2] + "\n")

    # pipe out ftp_partials to get_assembly_seqs_from_ftp.py
    # write to csv file?
    print(len(ftp_partials), ftp_partials)
