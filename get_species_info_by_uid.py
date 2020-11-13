# added to prevent error (comment out to see it)
from gevent import monkey as curious_george
curious_george.patch_all(thread=False, select=False)

from Bio import Entrez
from time import sleep
import re
import requests
import argparse


def get_protein_info(uniprot_ids):
    """
    Retrieves EMBL accession numbers and taxonomy ids for list of proteins. Creates a dict to map each protein's
    uid to its EMBL accession number and tax id.
    :param uniprot_ids: List of Uniprot IDs, e.g., ['P0AAJ3', 'A0NAQ1']
    :return: dictionary mapping each uid to its info
    """
    from bioservices import UniProt

    missing_embl = []
    missing_taxid = []
    orthos_map = {}
    u = UniProt()

    uniprot_records = list(map(lambda x: x.decode("utf-8"), u.retrieve(uniprot_ids, frmt='txt')))  # WSL CLI
    # uniprot_records = u.retrieve(uniprot_ids, frmt='txt')  # PyCharm

    embl_pattern = re.compile(r"DR\s+EMBL;.*?;\s+(.*?);")
    taxid_pattern = re.compile(r"OX\s+NCBI_TaxID=(\d+)")

    for i, record in enumerate(uniprot_records):
        embl_acc = get_match(embl_pattern, record, uniprot_ids, missing_embl, i)  # EMBL accession number for coding seq
        taxonomy_id = get_match(taxid_pattern, record, uniprot_ids, missing_taxid, i)  # tax_id of organism protein belongs to
        orthos_map[uniprot_ids[i]] = [embl_acc, taxonomy_id]  # map protein info to its uid

    if missing_embl:
        print('\n{} Protein(s) Missing EMBL Accession Number: '.format(len(missing_embl)) + ', '.join(missing_embl))

    if missing_taxid:
        print('\n{} Protein(s) Missing NCBI TaxID: '.format(len(missing_taxid)) + ', '.join(missing_taxid))

    return orthos_map


# rework/remove function?
def get_match(pattern, search_str, id_list, missing_list, index):
    """
    Find first match using RegEx; Allows for tracking of what IDs can't get the specific pattern match and placeholders
    to maintain indexing.
    :param pattern: compiled RegEx obj
    :param search_str: string to search
    :param id_list: list of identifiers
    :param missing_list: list to add identifier if pattern not found
    :param index: int; to show current spot in list of ids, for error handling only
    :return: Either the wanted pattern found (string), OR appends to list if error and returns placeholder None
    """

    match = re.findall(pattern, search_str)
    match = [x for x in match if x != '-']  # embl_acc are '-' if no Genomic DNA translation available

    if len(match) == 0:
        missing_list.append(id_list[index])
        # add placeholder value to keep indexing between refseq and uid 1:1, zipped later
        return None

    else:
        # grab first from list
        return match[0]


def get_assembly_id(tax_id):
    """
    Retrieves NCBI Assembly ID for a given NCBI Taxonomy ID. Email for Entrez must be set in main().
    :param tax_id: string; NCBI tax id
    :return: string; NCBI Assembly ID
    """
    search_term = f'txid{tax_id}'  # search term format used by NCBI Genome db
    handle1 = Entrez.esearch(db='genome', term=search_term)  # gets internal ID for genome in genome db
    record1 = Entrez.read(handle1)
    # if no info could be pulled for tax_id, don't want this protein... do try/except instead?
    if record1["Count"] == '0':
        return None
    else:
        org_id = record1['IdList'][0]    # grabs first id (more than 1 has not been encountered yet)

    handle2 = Entrez.esummary(db='genome', id=org_id)  # gives assembly of genome org
    record2 = Entrez.read(handle2)
    # if nothing returned from database (results in empty list)
    if not record2:
        return None
    else:
        assem_id = record2[0]["AssemblyID"]

    return assem_id


if __name__ == "__main__":

    # no optional args, but kept for clarity since there are a lot of arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType('r'), help="Fasta file containing protein UniProt IDs as header")
    parser.add_argument("outfile", type=argparse.FileType('w'), help="Name of fasta file to write coding seqs to")
    parser.add_argument('error_file', type=argparse.FileType('w'), help="Name of error log file")
    ftp = parser.add_argument_group('ftp', 'Arguments used to create csv file with partial ftp links')
    ftp.add_argument("ftp_out", type=argparse.FileType('w'), help="Filename for csv file with partial ftp links")
    ftp.add_argument('email', help='Email address required to use Entrez lookups')

    args = parser.parse_args()
    in_fh = args.infile
    fasta_fh = args.outfile
    error_fh = args.error_file   # can have duplicate uids logged, doesnt matter
    csv_fh = args.ftp_out
    Entrez.email = args.email

    # get list of uids
    uids = [x.strip().replace(">", "") for x in in_fh.readlines() if x.startswith(">")]

    ortho_mapping = get_protein_info(uids)

    # remove placeholder Nones
    tax_ids = [x[1] for x in ortho_mapping.values() if x[1] is not None]
    unique_tax_ids = list(set(tax_ids))  # need list for get_match() (set not subscriptable); if removed/reworked try just set

    # get partial ftp links to be used by get_assembly_seqs_from_ftp.py
    ftp_pattern = re.compile(r"<FtpPath_RefSeq>\S+(/genomes\S+GCF_\S+)<")
    ftp_partials = []
    missing_ftp = []
    for i, tax_id in enumerate(unique_tax_ids):
        assembly_id = get_assembly_id(tax_id)
        if assembly_id is None:
            missing_ftp.append(tax_id)
            continue

        handle = Entrez.esummary(db='assembly', id=assembly_id, report='full')
        string = handle.read().decode("utf-8")  # for WSL CLI
        # string = handle.read()    # Pycharm

        # gets ftp (partial) link for all cds in organism
        ftp_link = get_match(ftp_pattern, string, unique_tax_ids, missing_ftp, i)
        if ftp_link is None:
            sleep(0.4)  # prevent timeouts from ncbi
            continue    # already added to list by get_match()
        else:
            ftp_partials.append((ftp_link, tax_id))

        sleep(0.4)  # prevent timeouts from ncbi

    if missing_ftp:
        print('\n{} Organism(s) Missing FTP Link for Assembly: '.format(len(missing_ftp)) + ', '.join(missing_ftp))
        # error_fh.write('\n'.join(missing_ftp) + '\n')

    # replace tax_id with None for orthos that belong to orgs with unavailable assembly/ftp info (removed later)
    for uid, info in ortho_mapping.items():
        if info[1] in missing_ftp:
            error_fh.write(uid + '\n')
            info[1] = None

    # write to csv file
    for link, tax_id in ftp_partials:
        csv_fh.write(f'{link}, {tax_id}\n')
    csv_fh.close()

    # get CDS for gene and write all info to fasta file
    embl_base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
    for uid, info in ortho_mapping.items():
        # dont want to keep any ortho with incomplete info. Store uids in error file
        if None in info:
            error_fh.write(uid + '\n')    # add this to get_protein_info later?
            continue

        else:
            embl_acc = info[0]
            r = requests.get(embl_base_url + embl_acc)
            cds = ''.join(r.content.decode('utf-8').split("\n")[1:])
            # account for any http errors... will remove ortho from data for now; write uid info to std error?
            if 'status=' in cds or len(cds) == 0:
                print(f"{uid} removed from orthologs due to error in retrieving coding sequence.")
                error_fh.write(uid + '\n')
                continue

            else:
                fasta_fh.write(f">uid={uid};embl_acc={info[0]};tax_id={info[1]}\n")
                fasta_fh.write(cds + "\n")

    fasta_fh.close()
    error_fh.close()
